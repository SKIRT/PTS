#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.remoteframe Contains the RemoteFrame class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import tempfile
from itertools import count, izip

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.basics.remote import Remote, connected_remotes
from .frame import Frame
from .image import Image
from ...core.tools import parsing

# -----------------------------------------------------------------

def prepare_remote(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    # Check whether we are already connected to the specified remote host
    if host_id in connected_remotes and connected_remotes[host_id] is not None:
        remote = connected_remotes[host_id]
    else:

        # Debugging
        log.debug("Logging in to remote host ...")

        # Create a remote instance for the specified host ID
        remote = Remote()
        remote.setup(host_id)

    # Initiate python session
    if not remote.in_python_session: remote.start_python_session()

    # Import
    import_necessary_modules(remote)

    # Return the remote instance
    return remote

# -----------------------------------------------------------------

def import_necessary_modules(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    log.info("Importing necessary modules ...")

    # Import standard modules
    remote.import_python_package("tempfile")
    remote.import_python_package("urllib")
    remote.import_python_package("numpy", as_name="np")

    # Import the necessary PTS classes and modules
    remote.import_python_package("Frame", from_name="pts.magic.core.frame")
    remote.import_python_package("Image", from_name="pts.magic.core.image")
    remote.import_python_package("ConvolutionKernel", from_name="pts.magic.core.kernel")
    remote.import_python_package("CoordinateSystem", from_name="pts.magic.basics.coordinatesystem")
    remote.import_python_package("filesystem", from_name="pts.core.tools", as_name="fs")
    remote.import_python_package("archive", from_name="pts.core.tools")
    remote.import_python_package("parsing", from_name="pts.core.tools")

# -----------------------------------------------------------------

def get_first_missing_integer(integers):

    """
    This function ...
    :param integers: MUST BE SORTED !!!!
    :return:
    """

    if integers[0] != 0: return 0

    nums = (b for a, b in izip(integers, count(integers[0])) if a != b)
    return next(nums, integers[-1] + 1)

# -----------------------------------------------------------------

def get_frame_labels(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    variables = remote.python_variables()

    frame_labels = []

    for variable in variables:
        if variable.startswith("frame"):
            frame_labels.append(variable)

    return frame_labels

# -----------------------------------------------------------------

def get_frame_indices(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    return sorted([int(label.split("frame")[1]) for label in get_frame_labels(remote)])

# -----------------------------------------------------------------

def get_new_frame_label(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    current_indices = get_frame_indices(remote)
    return "frame" + str(get_first_missing_integer(current_indices))

# -----------------------------------------------------------------

def get_image_labels(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    variables = remote.python_variables()

    image_labels = []

    for variable in variables:
        if variable.startswith("image"):
            image_labels.append(variable)

    return image_labels

# -----------------------------------------------------------------

def get_image_indices(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    return sorted([int(label.split("image")[1]) for label in get_image_labels(remote)])

# -----------------------------------------------------------------

def get_new_image_label(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    current_indices = get_image_indices(remote)
    return "image" + str(get_first_missing_integer(current_indices))

# -----------------------------------------------------------------

#class MagicRemote(Remote):

    #"""
    #This class ...
    #"""

    #def __init__(self):

    #def setup(self):

# -----------------------------------------------------------------

class RemoteFrame(object):

    """
    This class ...
    """

    def __init__(self, label, remote):

        """
        This function ...
        :param label
        :param remote:
        """

        self.label = label
        self.remote = remote

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        output = self.remote.send_python_line("str(" + self.label + ")", output=True)
        return "".join(output)

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        output = self.remote.send_python_line("repr(" + self.label + ")", output=True)
        return "".join(output)

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Delete the Frame on the remote from the global namespace
        self.remote.remove_python_variable(self.label)

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = parsing.quantity(self.remote.get_simple_python_property(self.label, "fwhm"))
        return fwhm

    # -----------------------------------------------------------------

    @fwhm.setter
    def fwhm(self, fwhm):

        """
        This function ...
        :return:
        """

        self.remote.send_python_line(self.label + ".fwhm = parsing.quantity(" + str(fwhm) + ")")

    # -----------------------------------------------------------------

    @classmethod
    def from_url(cls, url, host_id, index=None, name=None, description=None, plane=None, hdulist_index=None,
                 no_filter=False, fwhm=None, add_meta=False):

        """
        This function ...
        :param url:
        :param host_id:
        :param index:
        :param name:
        :param description:
        :param plane:
        :param hdulist_index:
        :param no_filter:
        :param fwhm:
        :param add_meta:
        :return:
        """

        # Get remote instance
        remote = prepare_remote(host_id)

        # Inform the user
        log.info("Downloading file " + url + " ...")

        # Local path
        remote.send_python_line("temp_path = tempfile.gettempdir()")
        remote.send_python_line("filename = fs.name(url)")
        remote.send_python_line("local_path = fs.join(temp_path, filename)")

        # Download
        remote.send_python_line("urllib.urlretrieve(url, local_path)")

        # Get local path (on remote)
        local_path = remote.get_python_string("local_path")

        if local_path.endswith(".fits"): remote.send_python_line("fits_path = local_path")
        else:

            # Inform the user
            log.info("Decompressing kernel file ...")

            # Fits path
            remote.send_python_line("fits_path = fs.join(temp_path, fs.strip_extension(filename))")

            # Decompress the kernel FITS file
            remote.send_python_line("archive.decompress_file(local_path, fits_path)")

            # Remove the compressed file
            remote.send_python_line("fs.remove_file(local_path)")

        # Find label
        label = get_new_frame_label(remote)

        # Create RemoteFrame instance
        remoteframe = cls(label, remote)

        # Actually create the frame remotely
        remote.send_python_line(label + " = Frame.from_file(fits_path)")

        # Return the remoteframe instance
        return remoteframe

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, host_id, index=None, name=None, description=None, plane=None, hdulist_index=None,
                  no_filter=False, fwhm=None, add_meta=False):

        """
        This function ...
        :param path:
        :param host_id:
        :param index:
        :param name:
        :param description:
        :param plane:
        :param hdulist_index: if None, is automatically decided based on where the imageHDU is.
        :param no_filter:
        :param fwhm:
        :param add_meta:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Get the remote instance
        remote = prepare_remote(host_id)

        # Remote temp path
        remote.send_python_line("temp_path = tempfile.gettempdir()")
        remote_temp_path = remote.get_python_string("temp_path")

        # Get file name
        filename = fs.name(path)

        # Upload the frame file
        remote_frame_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the frame to " + remote_frame_path + " ...")
        remote.upload(path, remote_temp_path, compress=True, show_output=True)

        # Find label
        label = get_new_frame_label(remote)

        # Create RemoteFrame instance
        remoteframe = cls(label, remote)

        # Open the frame remotely
        log.info("Loading the frame on the remote host ...")

        # Actually create the frame remotely
        remote.send_python_line(label + " = Frame.from_file('" + remote_frame_path + "')")

        # Return the remoteframe instance
        return remoteframe

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Copying the remote frame ...")

        # Find new remote label
        label = get_new_frame_label(self.remote)

        # Copy the frame remotely
        self.remote.send_python_line(label + " = " + self.label + ".copy()")

        # Create new RemoteFrame instance
        newremoteframe = RemoteFrame(label, self.remote)

        # Return the new remoteframe instance
        return newremoteframe

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "xsize")

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "ysize")

    # -----------------------------------------------------------------

    @property
    def fwhm_pix(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "fwhm_pix")

    # -----------------------------------------------------------------

    def sum(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "sum()")

    # -----------------------------------------------------------------

    def normalize(self, to=1.0):

        """
        This function ...
        :param to:
        :return:
        """

        # Call the function on the remote
        self.remote.send_python_line(self.label + ".normalize(" + str(to) + ")")

    # -----------------------------------------------------------------

    def is_constant(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "is_constant()")

    # -----------------------------------------------------------------

    def convolved(self, *args, **kwargs):

        """
        This function ...
        :return:
        """

        new = self.copy()
        new.convolve(*args, **kwargs)
        return new

    # -----------------------------------------------------------------

    def convolve(self, kernel, allow_huge=False, fft=True):

        """
        This function ...
        :param kernel:
        :param allow_huge:
        :param fft:
        :return:
        """

        # Get the kernel FWHM
        kernel_fwhm = kernel.fwhm

        # Skip the calculation for a constant frame
        if self.is_constant():
            self.fwhm = kernel_fwhm
            return

        # SAVE KERNEL LOCALLY

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "kernel.fits")

        # Save the kernel locally
        kernel.save(local_path)

        # UPLOAD KERNEL

        # Remote temp path
        self.remote.send_python_line("temp_path = tempfile.gettempdir()")
        remote_temp_path = self.remote.get_python_string("temp_path")

        # Upload the kernel file
        remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
        log.info("Uploading the kernel to " + remote_kernel_path + " ...")
        self.remote.upload(local_path, remote_temp_path, compress=True, show_output=True)

        # Open the kernel remotely
        self.remote.send_python_line("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')")

        # Prepare the kernel if necessary
        #self.remote.send_python_line("if not kernel.prepared: kernel.prepare_for(" + self.label + ")")

        # Check where the NaNs are at
        #self.remote.send_python_line("nans_mask = np.isnan(" + self.label + "._data)")

        # Assert that the kernel is normalized
        #self.remote.send_python_line("assert kernel.normalized")

        # Do the convolution
        #self.remote. .... ETC


        # OR JUST: (does constant check again and creates copy for the prepared kernel, but... much more maintainable in this way)

        # Do the convolution on the remote frame
        self.remote.send_python_line(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ", fft=" + str(fft) + ")")

    # -----------------------------------------------------------------

    def rebinned(self, reference_wcs):

        """
        This function ...
        :param reference_wcs:
        :return:
        """

        new = self.copy()
        new.rebin(reference_wcs)
        return new

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=True, parallel=True):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :return:
        """

        # Upload the WCS
        self.remote.send_python_line("reference_wcs = CoordinateSystem('" + reference_wcs.to_header_string() + "')")

        # Create remote frame label for the footprint
        footprint_label = get_new_frame_label(self.remote)

        # Rebin, get the footprint
        self.remote.send_python_line(footprint_label + " = " + self.label + ".rebin(reference_wcs, exact=" + str(exact) + ", parallel=" + str(parallel) + ")")

        # Create a new remoteframe instance for the footprint
        remote_footprint_frame = RemoteFrame(footprint_label, self.remote)

        # Return the remoteframe footprint
        return remote_footprint_frame

    # -----------------------------------------------------------------

    def cropped(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        new = self.copy()
        new.crop(x_min, x_max, y_min, y_max)
        return new

    # -----------------------------------------------------------------

    def crop(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        # Crop remotely
        self.remote.send_python_line(self.label + ".crop(" + str(x_min) + ", " + str(x_max) + ", " + str(y_min) + ", " + str(y_max))

    # -----------------------------------------------------------------

    def padded(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        new = self.copy()
        new.pad(nx, ny)
        return new

    # -----------------------------------------------------------------

    def pad(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        # Pad remotely
        self.remote.send_python_line(self.label + ".pad(" + str(nx) + ", " + str(ny) + ")")

    # -----------------------------------------------------------------

    def unpad(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        self.remote.send_python_line(self.label + ".unpad(" + str(nx) + ", " + str(ny) + ")")

    # -----------------------------------------------------------------

    def downsampled(self, factor, order=3):

        """
        This function ...
        :return:
        """

        new = self.copy()
        new.downsample(factor, order)
        return new

    # -----------------------------------------------------------------

    def downsample(self, factor, order=3):

        """
        This function ...
        :param factor:
        :param order:
        :return:
        """

        self.remote.send_python_line(self.label + ".downsample(" + str(factor) + ", " + str(order) + ")")

    # -----------------------------------------------------------------

    def upsampled(self, factor, integers=False):

        """
        This function ...
        :param factor:
        :param integers:
        :return:
        """

        new = self.copy()
        new.upsample(factor, integers)
        return new

    # -----------------------------------------------------------------

    def upsample(self, factor, integers=False):

        """
        This function ...
        :param factor:
        :param integers:
        :return:
        """

        self.remote.send_python_line(self.label + ".upsample(" + str(factor) + ", " + str(integers) + ")")

    # -----------------------------------------------------------------

    def fill(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.remote.send_python_line(self.label + ".fill(" + str(value) + ")")

    # -----------------------------------------------------------------

    def replace_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.remote.send_python_line(self.label + ".replace_nans(" + str(value) + ")")

    # -----------------------------------------------------------------

    def replace_infs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.remote.send_python_line(self.label + ".replace_infs(" + str(value) + ")")

    # -----------------------------------------------------------------

    @classmethod
    def from_local(cls, frame, host_id):

        """
        This function ...
        :param frame:
        :param host_id:
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "newframe_fromlocal.fits")

        # Save the frame locally
        frame.save(local_path)

        # Create the remoteframe from the locally saved frame
        remoteframe = cls.from_file(local_path, host_id)

        # Remove the local file
        fs.remove_file(local_path)

        # Return the new remoteframe
        return remoteframe

    # -----------------------------------------------------------------

    def to_local(self):

        """
        This function ...
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, self.label + ".fits")

        # Save the remote frame locally to the temp path
        self.save(local_path)

        # Create the frame
        frame = Frame.from_file(local_path)

        # Remove the local temporary file
        fs.remove_file(local_path)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        """

        # Inform the user
        log.info("Saving the remote frame to '" + path + "' locally ...")

        # Determine filename and local directory
        filename = fs.name(path)
        local_directory = fs.directory_of(path)

        # Determine remote temp path
        self.remote.send_python_line("temp_path = tempfile.gettempdir()")

        # Debugging
        log.debug("Saving the frame remotely ...")

        # Determine path to save the frame remotely first
        self.remote.send_python_line("remote_path = fs.join(temp_path, " + filename + ")")
        remote_path = self.remote.get_python_string("remote_path")

        # Debugging
        log.debug("Downloading the frame ...")

        # Save the frame remotely
        self.remote.send_python_line(self.label + ".save(remote_path)")

        # Download
        self.remote.download(remote_path, local_directory)

        # Remove the remote file
        self.remote.remove_file(remote_path)

# -----------------------------------------------------------------

class RemoteImage(object):

    """
    This class ...
    """

    def __init__(self, label, remote):

        """
        The constructor ...
        :param label:
        :param remote:
        """

        self.label = label
        self.remote = remote

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        output = self.remote.send_python_line("str(" + self.label + ")", output=True)
        return "".join(output)

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        output = self.remote.send_python_line("repr(" + self.label + ")", output=True)
        return "".join(output)

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Delete the Image on the remote from the global namespace
        self.remote.remove_python_variable(self.label)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, host_id):

        """
        This function ...
        :param path:
        :param host_id:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Get the remote instance
        remote = prepare_remote(host_id)

        # Remote temp path
        remote.send_python_line("temp_path = tempfile.gettempdir()")
        remote_temp_path = remote.get_python_string("temp_path")

        # Get file name
        filename = fs.name(path)

        # Upload the image file
        remote_image_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the image to " + remote_image_path + " ...")
        remote.upload(path, remote_temp_path, compress=True, show_output=True)

        # Find label
        label = get_new_image_label(remote)

        # Create RemoteImage instance
        remoteimage = cls(label, remote)

        # Open the frame remotely
        log.info("Loading the image on the remote host ...")

        # Actually create the frame remotely
        remote.send_python_line(label + " = Image.from_file('" + remote_image_path + "')")

        # Return the remoteimage instance
        return remoteimage

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the remote image to '" + path + "' locally ...")

        # Determine filename and local directory
        filename = fs.name(path)
        local_directory = fs.directory_of(path)

        # Determine remote temp path
        self.remote.send_python_line("temp_path = tempfile.gettempdir()")

        # Debugging
        log.debug("Saving the image remotely ...")

        # Determine path to save the frame remotely first
        self.remote.send_python_line("remote_path = fs.join(temp_path, " + filename + ")")
        remote_path = self.remote.get_python_string("remote_path")

        # Debugging
        log.debug("Downloading the image ...")

        # Save the image remotely
        self.remote.send_python_line(self.label + ".save(remote_path)")

        # Download
        self.remote.download(remote_path, local_directory)

        # Remove the remote file
        self.remote.remove_file(remote_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_local(cls, image, host_id):

        """
        This function ...
        :param image:
        :param host_id:
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "newimage_fromlocal.fits")

        # Save the image locally
        image.save(local_path)

        # Create the remoteimage from the locally saved image
        remoteimage = cls.from_file(local_path, host_id)

        # Remove the local file
        fs.remove_file(local_path)

        # Return the new remoteimage
        return remoteimage

    # -----------------------------------------------------------------

    def to_local(self):

        """
        This function ...
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, self.label + ".fits")

        # Save the remote image locally to the temp path
        self.save(local_path)

        # Create the image
        image = Image.from_file(local_path)

        # Remove the local temporary file
        fs.remove_file(local_path)

        # Return the image
        return image

    # -----------------------------------------------------------------

    def convolve(self, kernel, allow_huge=False):

        """
        This function ...
        :param kernel:
        :param allow_huge:
        :return:
        """

        # SAVE KERNEL LOCALLY

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "kernel.fits")

        # Save the kernel locally
        kernel.save(local_path)

        # UPLOAD KERNEL

        # Remote temp path
        self.remote.send_python_line("temp_path = tempfile.gettempdir()")
        remote_temp_path = self.remote.get_python_string("temp_path")

        # Upload the kernel file
        remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
        log.info("Uploading the kernel to " + remote_kernel_path + " ...")
        self.remote.upload(local_path, remote_temp_path, compress=True, show_output=True)

        # Open the kernel remotely
        self.remote.send_python_line("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')")

        # Convolve the image remotely
        self.remote.send_python_line(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ")")

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs):

        """
        This function ...
        :return:
        """

        # Upload the WCS
        self.remote.send_python_line("reference_wcs = CoordinateSystem('" + reference_wcs.to_header_string() + "')")

        # Rebin remotely
        self.remote.send_python_line(self.label + ".rebin(reference_wcs)")

# -----------------------------------------------------------------
