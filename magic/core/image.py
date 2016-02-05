#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.image Contains the Image class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np
import matplotlib.pyplot as plt
import copy

# Import astronomical modules
import aplpy
from astropy.io import fits
from astropy import units as u

# Import the relevant AstroMagic classes and modules
from ..basics import Layers, Region, Mask, CoordinateSystem
from .frame import Frame
from ..tools import headers, fitting, plotting, statistics, catalogs, transformations
from ...core.tools.logging import log

# -----------------------------------------------------------------

class Image(object):

    """
    This class ...
    """

    def __init__(self, filepath=None, always_call_first_primary=True):

        """
        The constructor ...
        :param filepath:
        :return:
        """

        # Initialize a set of layers to represent image frames, masks and regions
        self.frames = Layers()
        self.masks = Layers()
        self.regions = Layers()

        # The image name and path
        self.name = None
        self.path = None

        # The original image header
        self.original_header = None

        # The dictionary containing meta information
        self.metadata = dict()

        if filepath is not None:

            # Check if the specified file exists, otherwise exit with an error
            if not os.path.isfile(filepath): raise IOError("No such file: " + filepath)

            # Set the name of the image
            self.name = os.path.splitext(os.path.basename(filepath))[0]

            # Read in the image
            self.load_frames(filepath, always_call_first_primary=always_call_first_primary)
            self.path = filepath

    # -----------------------------------------------------------------

    def select_all(self):

        """
        This function ...
        :return:
        """

        # Select all frames, regions and masks
        self.frames.select_all()
        self.regions.select_all()
        self.masks.select_all()

    # -----------------------------------------------------------------

    def deselect_all(self):

        """
        This function ...
        :return:
        """

        # Deselect all frames, regions and masks
        self.frames.deselect_all()
        self.regions.deselect_all()
        self.masks.deselect_all()

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None
        return self.frames.primary.shape

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None
        return self.frames.primary.xsize

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None
        return self.frames.primary.ysize

    # -----------------------------------------------------------------

    @property
    def filter(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None

        # Return the filter of the primary frame
        return self.frames.primary.filter

    # -----------------------------------------------------------------

    @property
    def wavelength(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None

        # Return the wavelength of the primary frame
        return self.frames.primary.wavelength

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None

        # Return the unit of the primary frame
        return self.frames.primary.unit

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None

        # Return the pixelscale of the primary frame
        return self.frames.primary.pixelscale

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None

        # Return the FWHM of the primary frame
        return self.frames.primary.fwhm

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        if "primary" not in self.frames: return None

        # Return the wcs of the primary frame
        return self.frames.primary.wcs

    # -----------------------------------------------------------------

    def save(self, path=None, add_metadata=False, origin=None):

        """
        This function exports the currently selected frame(s) as a datacube into FITS file
        :param path:
        :return:
        """

        if path is None: path = self.path

        # Create an array to contain the data cube
        datacube = []

        plane_index = 0

        header = None

        # Export all frames to the specified file
        for frame_name in self.frames:

            # Inform the user that this frame will be saved to the image file
            log.info("Exporting the " + frame_name + " frame to " + path)

            if header is None: header = self.frames[frame_name].header

            # Check if the coordinate system of this frame matches that of the other frames
            #if header != self.frames[frame_name].header: raise ValueError("The WCS of the different frames does not match")

            # Add this frame to the data cube, if its coordinates match those of the primary frame
            datacube.append(self.frames[frame_name])
            
            # Add the name of the frame to the header
            header["PLANE" + str(plane_index)] = frame_name + " [frame]"

            # Increment the plane index
            plane_index += 1

        # Export all masks to the specified file
        for mask_name in self.masks:

            # Inform the user that this mask will be saved to the image file
            log.info("Exporting the " + mask_name + " mask to " + path)

            # Add this mask to the data cube
            datacube.append(self.masks[mask_name].astype(int))

            # Add the name of the mask to the header
            header["PLANE" + str(plane_index)] = mask_name + " [mask]"

            # Increment the plane index
            plane_index += 1

        # Add the meta information to the header
        if add_metadata:
            for key in self.metadata:
                try: header[key] = self.metadata[key]
                except ValueError: pass # Some values in the header gives errors in Astropy when adding them again to this new header ... (e.g. ValueError: Illegal value: = 'created by T.H. Jarrett'.)

        # Set plane information
        if plane_index > 1:
            header["NAXIS"] = 3
            header["NAXIS3"] = plane_index

        # Set the unit
        header["BUNIT"] = str(self.unit)

        # Add origin description
        if origin is not None: header["ORIGIN"] = origin
        else: header["ORIGIN"] = "Image class of PTS package"

        # Create the HDU from the data array and the header
        hdu = fits.PrimaryHDU(np.array(datacube), header)

        # Write the HDU to a FITS file
        hdu.writeto(path, clobber=True)

        # Update the path
        self.path = path

        # Inform the user that the file has been created
        log.info("File " + path + " created")

    # -----------------------------------------------------------------

    def import_region(self, path, name, overwrite=False):

        """
        This function imports a new region from a DS9 region file
        :param path:
        :param name:
        :return:
        """

        # Create an Region object from the regions file
        region = Region.from_file(path)

        # Add the region to the set of regions
        self.add_region(region, name, overwrite)

    # -----------------------------------------------------------------

    def export_region(self, path):

        """
        This function exports the currently selected region(s) to a DS9 region file
        :param path:
        :return:
        """

        # Find the active region
        region = self.regions.get_selected(require_single=True)

        # Inform the user
        log.info("Creating " + path + " from the " + region + " region")

        # Write the region file
        self.regions[region].region.write(path)

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

        # Deselect all frames, regions and masks of this image
        self.deselect_all()

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

    def set_unit(self, unit):

        """
        This function ...
        """

        # Loop over all currently selected frames
        for frame_name in self.frames.get_selected():

            # Inform the user
            log.info("Setting the unit of the " + frame_name + " frame to " + str(unit))

            # Set the unit for this frame
            self.frames[frame_name].set_unit(unit)

    # -----------------------------------------------------------------

    def set_fwhm(self, fwhm):

        """
        This function ...
        :param fwhm:
        :return:
        """

        # Loop over all currently selected frames
        for frame_name in self.frames.get_selected():

            # Inform the user
            log.info("Setting the FWHM of the " + frame_name + " frame to " + str(fwhm))

            # Set the unit for this frame
            self.frames[frame_name].set_fwhm(fwhm)

    # -----------------------------------------------------------------

    def convert_to(self, unit):

        """
        This function ...
        """

        # Loop over all currently selected frames
        for frame_name in self.frames.get_selected():

            # Inform the user
            log.info("Converting the unit of the " + frame_name + " frame to " + str(unit))

            # Set the unit for this frame
            self.frames[frame_name].set_unit(unit)

    # -----------------------------------------------------------------

    def convolve(self, kernel):

        """
        This function ...
        """

        # Loop over all currently selected frames
        for frame_name in self.frames:

            # Inform the user
            log.info("Convolving the " + frame_name + " frame ...")

            # Convolve this frame
            self.frames[frame_name] = self.frames[frame_name].convolved(kernel)

    # -----------------------------------------------------------------

    def rebin(self, reference):

        """
        This function ...
        """

        # Loop over all currently selected frames
        for frame_name in self.frames:

            # Inform the user
            log.info("Rebinning the " + frame_name + " frame ...")

            # Rebin this frame (the reference wcs is automatically set in the new frame)
            self.frames[frame_name] = self.frames[frame_name].rebinned(reference)

        # Loop over the masks
        for mask_name in self.masks:

            # Inform the user
            log.info("Rebinning the " + mask_name + " mask ...")

            # Rebin this mask
            data = transformations.align_and_rebin(self, self.wcs.to_header(), reference.header)

            # Return the rebinned mask
            self.masks[mask_name] = Mask(data, self.masks[mask_name].selected, self.masks[mask_name].description)

    # -----------------------------------------------------------------

    def crop(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        """

        # Loop over all currently selected frames
        for frame_name in self.frames.get_selected():

            # Inform the user
            log.info("Cropping the " + frame_name + " frame")

            # Rebin this frame
            self.frames[frame_name] = self.frames[frame_name].crop(x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    def plot(self, path=None, color=True, grid=False, blacknan=False, publication=False):

        """
        This function shows a plot of the currently selected frame, combined with the active regions and masks
        :param path:
        :param color:
        :param grid:
        :param blacknan:
        :param publication:
        :return:
        """

        # Get the currently active frame
        frame = self.frames.get_selected()[0]

        # Create a total mask of the currently active masks
        total_mask = self.combine_masks(return_mask=True)

        # Mask the frame with nans
        maskedimage = np.ma.array(self.frames[frame], mask = total_mask)
        image_with_nans =  maskedimage.filled(np.NaN)

        # Create a HDU from this frame with the image header
        hdu = fits.PrimaryHDU(image_with_nans, self.original_header)

        if path is None:

            # Create a figure canvas
            figure = plt.figure(figsize=(12, 12))

            # Create a figure from this frame
            plot = aplpy.FITSFigure(hdu, figure=figure)

        else:

            # Create a figure from this frame
            plot = aplpy.FITSFigure(hdu)

        if color:

            # Plot in color scale
            plot.show_colorscale()

        else:

            # Plot in gray scale
            plot.show_grayscale()

        # Add a color bar
        plot.add_colorbar()

        if blacknan:

            # Set the nan color to black
            plot.set_nan_color('black')

        if grid:

            # Add a grid
            plot.add_grid()

        # If requested, use the 'publication' theme
        if publication: plot.set_theme('publication')

        # Add the regions
        for region in self.regions.get_selected():

            # Get the shape list
            shapes = self.regions[region].region.as_imagecoord(self.original_header)

            # Add these shapes to the plot
            plot.show_regions(shapes)

        if path is None:

            #plt.draw()
            #plt.close('all') # redundant
            #plt.show(block=False)
            plt.show()

        else: plot.save(path)

    # -----------------------------------------------------------------

    def delete_frames(self):

        """
        This function removes the currently selected frame(s)
        :return:
        """

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Deleting the " + frame_name + " frame")

            # Remove this frame from the frames dictionary
            del self.frames[frame_name]

    # -----------------------------------------------------------------

    def copy_frames(self):

        """
        This function ...
        :return:
        """

        # For each selected frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Copying the " + frame_name + " frame as another frame")

            # Copy the data and add it as a new frame
            data_copy = copy.deepcopy(self.frames[frame_name])
            #coordinates = self.frames[frame_name].coordinates

            data_copy.description = "Copy of the "+frame_name+" frame"

            self.add_frame(data_copy, frame_name+"_copy")

    # -----------------------------------------------------------------

    def delete_regions(self):

        """
        This function removes the currently selected region(s)
        :return:
        """

        # For each active region
        for region_name in self.regions.get_selected(allow_none=False):

            # Inform the user
            log.info("Deleting the " + region_name + " region")

            # Remove this region from the regions dictionary
            del self.regions[region_name]

    # -----------------------------------------------------------------

    def delete_masks(self):

        """
        This function removes the currently selected mask(s)
        :return:
        """

        # For each active mask
        for mask_name in self.masks.get_selected(allow_none=False):

            # Inform the user
            log.info("Deleting the " + mask_name + " mask")

            # Remove this mask from the masks dictionary
            del self.masks[mask_name]

    # -----------------------------------------------------------------

    def apply_masks(self, fill=0.0):

        """
        This function ...
        :param fill:
        :return:
        """

        # Loop over all selected frames
        for frame_name in self.frames.get_selected(allow_none=False):

            # Loop over all selected masks
            for mask_name in self.masks.get_selected(allow_none=False):

                # Inform the user
                log.info("Applying the " + mask_name + " mask to the " + frame_name + " frame")

                # Apply the mask
                self.masks[mask_name].apply(self.frames[frame_name], fill)

    # -----------------------------------------------------------------

    def combine_masks(self, name=None, allow_none=True, return_mask=False):

        """
        This function ...
        :param name:
        :param allow_none:
        :return:
        """

        # TODO: rewrite this function
        pass

    # -----------------------------------------------------------------

    def __imul__(self, factor):

        """
        This function ...
        :param other:
        :return:
        """

        # Loop over all currently selected frames
        for frame_name in self.frames.get_selected():

            # Inform the user
            log.info("Multiplying the " + frame_name + " frame by a factor of " + str(factor))

            # Multiply the frame by the given factor
            self.frames[frame_name] *= factor

        # Return a reference to this instance
        return self

    # -----------------------------------------------------------------

    def __idiv__(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        # Loop over all currently selected frames
        for frame_name in self.frames.get_selected():

            # Inform the user
            log.info("Dividing the " + frame_name + " frame by a factor of " + str(factor))

            # Divide the frame by the given factor
            self.frames[frame_name] /= factor

        # Return a reference to this instance
        return self

    # -----------------------------------------------------------------

    def create_mask(self, return_mask=False):

        """
        This function creates a mask from the currently selected region(s)
        :return:
        """

        # TODO: rewrite this function
        pass

    # -----------------------------------------------------------------

    def get_galactic_extinction(self, galaxy_name):

        """
        This function ...
        """

        return catalogs.fetch_galactic_extinction(galaxy_name, self.filter)

    # -----------------------------------------------------------------

    def rename_region(self, name):

        """
        This function renames a region
        :param name:
        :return:
        """

        # Get the name of the currently selected region
        region_name = self.regions.get_selected(require_single=True)

        # Remove the region of the dictionary of regions and re-add it under a different key
        self.regions[name] = self.regions.pop(region_name)

    # -----------------------------------------------------------------

    def rename_frame(self, name):

        """
        This function renames a frame
        :param name:
        :return:
        """

        # Get the name of the currently selected frame
        frame_name = self.frames.get_selected(require_single=True)

        # Remove the frame from the dictionary of frames and re-add it under a different key
        self.frames[name] = self.frames.pop(frame_name)

    # -----------------------------------------------------------------

    def rename_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the name of the currently selected mask
        mask_name = self.masks.get_selected(require_single=True)

        # Remove the mask from the dictionary of masks and re-add it under a different key
        self.masks[name] = self.masks.pop(mask_name)

    # -----------------------------------------------------------------

    def fit_polynomial(self, plot=False, degree=3, sigma_clipping=True):

        """
        This function fits a polynomial function to each of the currently active frames
        :param plot:
        :param upsample_factor:
        :return:
        """

        # Get the currently active mask
        total_mask = self.combine_masks(return_mask=True)

        # For each currently active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Fitting a polynomial function to the " + frame_name + " frame")

            if sigma_clipping: new_mask = statistics.sigma_clip_mask(self.frames[frame_name], mask=total_mask)
            else: new_mask = total_mask

            # Fit the model
            polynomial = fitting.fit_polynomial(self.frames[frame_name], mask=new_mask, degree=degree)

            # Plot the difference between the data and the model, if requested
            if plot: plotting.plot_difference_model(self.frames[frame_name], polynomial)

            # Evaluate the model
            evaluated = fitting.evaluate_model(polynomial, 0, self.frames[frame_name].xsize, 0, self.frames[frame_name].ysize)

            # Add the evaluated model as a new frame
            description = "A polynomial fit to the " + frame_name + " primary frame"
            self.add_frame(Frame(evaluated, self.frames[frame_name].coordinates, description), frame_name+"_polynomial")

    # -----------------------------------------------------------------

    def load_frames(self, filename, index=None, name=None, description=None, always_call_first_primary=True, rebin_to_wcs=False):

        """
        This function ...
        :param filename:
        :param index:
        :param name:
        :param description:
        :param always_call_first_primary:
        :param rebin_to_wcs:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file: " + filename)

        # Open the HDU list for the FITS file
        hdulist = fits.open(filename)

        # Get the primary HDU
        hdu = hdulist[0]

        # Get the image header
        header = hdu.header

        # Get a copy of the original header
        original_header = header.copy()

        # Remove references to the third axis
        header["NAXIS"] = 2
        if "NAXIS3" in header: del header["NAXIS3"]
        for key in header:
            if "PLANE" in key: del header[key]

        # Obtain the world coordinate system
        wcs = CoordinateSystem(header)

        # Pixelscale: WILL BE CHECKED IN THE FRAME CONSTRUCTOR

        # Obtain the filter for this image
        filter = headers.get_filter(self.name, original_header)

        # Inform the user on the filter
        if filter is not None: log.info("The filter for this image is " + filter.filterID())
        else: log.warning("Could not determine the filter for this image")

        # Obtain the units of this image
        unit = headers.get_unit(original_header)

        # Get the magnitude zero-point
        zero_point = headers.get_zero_point(header)

        # Check whether the image is sky-subtracted
        sky_subtracted = headers.is_sky_subtracted(original_header)

        # Check whether multiple planes are present in the FITS image
        nframes = headers.get_number_of_frames(original_header)
        if nframes > 1:

            # For each frame
            for i in range(nframes):

                # If only a frame with specific index needs to be imported, skip this frame if it does not correspond
                if index is not None and i != index: continue

                if index is None:

                    # Get the name of this frame, but the first frame always gets the name 'primary' unless the
                    # 'always_call_first_primary' flag is disabled
                    if i == 0 and always_call_first_primary:
                        description = "the primary signal map"
                        name = "primary"
                    else:

                        description = headers.get_frame_description(original_header, i)

                        if description is None:
                            description = ""
                            name = "frame"+str(i)
                        else: name = headers.get_frame_name(description)

                # The sky-subtracted flag should only be set for the primary frame
                subtracted = sky_subtracted if i == 0 else False

                # Check the shape of this new frame
                if self.shape is not None:

                    if hdu.data[i].shape != self.shape:

                        if rebin_to_wcs:

                            # Inform the user
                            log.warning("Rebinning the " + name + " frame (plane " + str(i) + ") of " + filename + " to match the shape of this image")

                            # Check if the unit is a surface brightness unit
                            if unit != u.Unit("MJy/sr"): raise ValueError("Cannot rebin since unit " + str(unit) + " is not recognized as a surface brightness unit")

                            # Change the data and the WCS
                            hdu.data[i] = transformations.align_and_rebin(hdu.data[i], header, self.wcs.to_header())
                            wcs = self.wcs

                        else: raise ValueError("The shape of the " + name + " frame (plane " + str(i) + ") of " + filename + " does not match the shape of this image")

                # Add this frame to the frames dictionary
                self.add_frame(Frame(hdu.data[i], wcs, description, False, unit, name, filter, subtracted, zero_point), name)

                # Select the frame if it is the first one (the primary frame)
                if i == 0: self.frames[name].select()

        else:

            # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
            if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

            if name is None: name = "primary"
            if description is None: description = "the primary signal map"

            # Check the shape of this new frame
            if self.shape is not None:

                if hdu.data.shape != self.shape:

                    if rebin_to_wcs:

                        # Inform the user
                        log.warning("Rebinning the " + name + " frame (plane 0) of " + filename + " to match the shape of this image")

                        # Check if the unit is a surface brightness unit
                        if unit != u.Unit("MJy/sr"): raise ValueError("Cannot rebin since unit " + str(unit) + " is not recognized as a surface brightness unit")

                        # Change the data and the WCS
                        hdu.data = transformations.align_and_rebin(hdu.data, header, self.wcs.to_header())
                        wcs = self.wcs

                    else: raise ValueError("The shape of the " + name + " frame (plane 0) of " + filename + " does not match the shape of this image")

            # Add the primary image frame
            self.add_frame(Frame(hdu.data, wcs, description, False, unit, name, filter, sky_subtracted, zero_point), name)

            # Select the frame
            self.frames[name].select()

        # Add meta information
        for key in original_header:
            self.metadata[key.lower()] = original_header[key]

        # Let the original header remain accessible
        self.original_header = original_header

        # Close the FITS file
        hdulist.close()

    # -----------------------------------------------------------------

    def add_frame(self, frame, name, overwrite=False):

        """
        This function ...
        :param frame:
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adding '" + name + "' to the set of frames ...")

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

    def add_region(self, region, name, overwrite=False):

        """
        This function ...
        :param region:
        :param name:
        :param overwrite:
        :return:
        """

        # Inform the user
        log.info("Adding '" + name + "' to the set of regions ...")

        # Check whether a region with this name already exists
        if name in self.regions and not overwrite: raise RuntimeError("A region with this name already exists")

        # Add the region to the set of regions
        self.regions[name] = region

    # -----------------------------------------------------------------

    def add_mask(self, mask, name, overwrite=False):

        """
        This function ...
        :param mask:
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adding '" + name + "' to the set of masks ...")

        # Check whether a mask with this name already exists
        if name in self.masks and not overwrite: raise RuntimeError("A mask with this name already exists")

        # Check if the shape matches the shape of this image
        if mask.shape != self.shape: raise ValueError("Mask does not have the correct shape for this image")

        # Add the mask to the set of masks
        self.masks[name] = mask

# -----------------------------------------------------------------
