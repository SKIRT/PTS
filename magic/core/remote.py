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
from .frame import Frame # IMPORTANT THAT THESE ARE IMPORTED !!
from .image import Image # IMPORTANT THAT THESE ARE IMPORTED !!
from .datacube import DataCube # IMPORTANT THAT THESE ARE IMPORTED !!
from ...core.tools import parsing

# -----------------------------------------------------------------

prepared = [] # the list of host IDs for which the remote has been prepared

# -----------------------------------------------------------------

def prepare_remote(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    # If this host has already been prepared, just return the prepared remote
    if host_id in prepared: return connected_remotes[host_id]

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

    # Add the host ID to the 'prepared' list
    prepared.append(host_id)

    # Return the remote instance
    return remote

# -----------------------------------------------------------------

def import_necessary_modules(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    # Inform the user
    log.info("Importing necessary modules ...")

    # Import standard modules
    #remote.import_python_package("tempfile")  ## doesn't work: we seem to have no permissions in this directory on nancy
    remote.import_python_package("urllib")
    remote.import_python_package("numpy", as_name="np")

    # Import the necessary PTS classes and modules
    remote.import_python_package("Frame", from_name="pts.magic.core.frame")
    remote.import_python_package("Image", from_name="pts.magic.core.image")
    remote.import_python_package("DataCube", from_name="pts.magic.core.datacube")
    remote.import_python_package("ConvolutionKernel", from_name="pts.magic.core.kernel")
    remote.import_python_package("CoordinateSystem", from_name="pts.magic.basics.coordinatesystem")
    remote.import_python_package("archive", from_name="pts.core.tools")
    remote.import_python_package("parsing", from_name="pts.core.tools")
    remote.import_python_package("Filter", from_name="pts.core.basics.filter")

# -----------------------------------------------------------------

def get_first_missing_integer(integers):

    """
    This function ...
    :param integers: MUST BE SORTED !!!!
    :return:
    """

    if len(integers) == 0: return 0

    if integers[0] != 0: return 0

    nums = (b for a, b in izip(integers, count(integers[0])) if a != b)
    return next(nums, integers[-1] + 1)

# -----------------------------------------------------------------

def get_labels(classname, remote):

    """
    This function ...
    :param classname:
    :param remote:
    :return:
    """

    variables = remote.python_variables()

    labels = []

    for variable in variables:
        if variable.startswith(classname.lower()):
            labels.append(variable)

    return labels

# -----------------------------------------------------------------

def get_indices(classname, remote):

    """
    This function ...
    :param classname:
    :param remote:
    :return:
    """

    return sorted([int(label.split(classname.lower())[1]) for label in get_labels(classname, remote)])

# -----------------------------------------------------------------

def get_new_label(classname, remote):

    """
    This function ...
    :param classname:
    :param remote:
    :return:
    """

    current_indices = get_indices(classname, remote)
    return classname.lower() + str(get_first_missing_integer(current_indices))

# -----------------------------------------------------------------

class RemoteFrame(object):

    """
    This class ...
    """

    @classmethod
    def local_class(cls):

        """
        This function ...
        :return:
        """

        classname = cls.local_classname
        return globals()[classname]

    # -----------------------------------------------------------------

    @classmethod
    def local_classname(cls):

        """
        This function ...
        :return:
        """

        return cls.__name__.split("Remote")[1]

    # -----------------------------------------------------------------

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

        return self.remote.get_simple_python_property(self.label, "__str__()")

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "__repr__()")

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
                 no_filter=False, fwhm=None, add_meta=True):

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
        remote_temp_path = remote.session_temp_directory
        remote.send_python_line("filename = fs.name(url)")
        remote.send_python_line("local_path = fs.join('" + remote_temp_path + "', filename)")

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
        label = get_new_label(cls.local_classname(), remote)

        # Create RemoteFrame instance
        remoteframe = cls(label, remote)

        # Actually create the frame remotely
        remote.send_python_line(label + " = " + cls.local_classname() + ".from_file(fits_path)")

        # Return the remoteframe instance
        return remoteframe

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, host_id, index=None, name=None, description=None, plane=None, hdulist_index=None,
                  no_filter=False, fwhm=None, add_meta=True):

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
        remote_temp_path = remote.session_temp_directory

        # Get file name
        filename = fs.name(path)

        # Upload the frame file
        remote_frame_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the frame to " + remote_frame_path + " ...")
        remote.upload(path, remote_temp_path, compress=True, show_output=True)

        # Find label
        label = get_new_label(cls.local_classname(), remote)

        # Create RemoteFrame instance
        remoteframe = cls(label, remote)

        # Open the frame remotely
        log.info("Loading the frame on the remote host ...")

        # Actually create the frame remotely
        remote.send_python_line(label + " = " + cls.local_classname() + ".from_file('" + remote_frame_path + "')")

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
        label = get_new_label(self.local_classname(), self.remote)

        # Copy the frame remotely
        self.remote.send_python_line(label + " = " + self.label + ".copy()")

        # Create new RemoteFrame instance
        newremoteframe = self.__class__(label, self.remote)

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
        remote_temp_path = self.remote.session_temp_directory

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
        self.remote.send_python_line(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ", fft=" + str(fft) + ")", show_output=True, timeout=None)

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
        footprint_label = get_new_label("Frame", self.remote)

        # Rebin, get the footprint
        self.remote.send_python_line(footprint_label + " = " + self.label + ".rebin(reference_wcs, exact=" + str(exact) + ", parallel=" + str(parallel) + ")", timeout=None)

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

        # Create the local frame
        frame = self.local_class().from_file(local_path)

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
        remote_temp_path = self.remote.session_temp_directory

        # Remote file path
        remote_file_path = fs.join(remote_temp_path, filename)

        # Debugging
        log.debug("Saving the frame remotely ...")

        # Save the frame remotely
        self.remote.send_python_line(self.label + ".save('" + remote_file_path + "')")

        # Debugging
        log.debug("Downloading the frame ...")

        # Download
        self.remote.download(remote_file_path, local_directory)

        # Remove the remote file
        self.remote.remove_file(remote_file_path)

# -----------------------------------------------------------------

class RemoteImage(object):

    """
    This class ...
    """

    @classmethod
    def local_class(cls):

        """
        This function ...
        :return:
        """

        classname = cls.local_classname
        return globals()[classname]

    # -----------------------------------------------------------------

    @classmethod
    def local_classname(cls):
        """
        This function ...
        :return:
        """

        return cls.__name__.split("Remote")[1]

    # -----------------------------------------------------------------

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

        return self.remote.get_simple_python_property(self.label, "__str__()")

        # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "__repr__()")

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
        remote_temp_path = remote.session_temp_directory

        # Get file name
        filename = fs.name(path)

        # Upload the image file
        remote_image_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the image to '" + remote_image_path + "' ...")
        remote.upload(path, remote_temp_path, compress=True, show_output=True)

        # Find label
        label = get_new_label(cls.local_classname(), remote)

        # Create RemoteImage instance
        remoteimage = cls(label, remote)

        # Open the frame remotely
        log.info("Loading the image on the remote host ...")

        # Actually create the frame remotely
        remote.send_python_line(label + " = " + cls.local_classname() + ".from_file('" + remote_image_path + "')")

        # Return the remoteimage instance
        return remoteimage

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Copying the remote image ...")

        # Find new remote label
        label = get_new_label(self.local_classname(), self.remote)

        # Copy the frame remotely
        self.remote.send_python_line(label + " = " + self.label + ".copy()")

        # Create new RemoteImage instance
        newremoteimage = self.__class__(label, self.remote)

        # Return the new remoteimage instance
        return newremoteimage

    # -----------------------------------------------------------------

    @property
    def primary(self):

        """
        This function ...
        :return:
        """

        # Assign a remote variable to the primary frame of this remote image
        label = get_new_label("Frame", self.remote)

        # Reference the primary frame to the new variable
        self.remote.send_python_line(label + " = " + self.label + ".primary")

        # Create a new RemoteFrame instance
        newremoteframe = RemoteFrame(label, self.remote) # We assume the frame of the remote image is not of a type derived from Frame

        # Return the new remoteframe instance
        return newremoteframe

    # -----------------------------------------------------------------

    @property
    def has_frames(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "has_frames")

    # -----------------------------------------------------------------

    @property
    def nframes(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "nframes")

    # -----------------------------------------------------------------

    @property
    def nmasks(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "nmasks")

    # -----------------------------------------------------------------

    @property
    def nregions(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "nregions")

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return self.remote.get_simple_python_property(self.label, "shape")

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
        remote_temp_path = self.remote.session_temp_directory

        # Determine remote path for the image
        remote_image_path = fs.join(remote_temp_path, filename)

        # Debugging
        log.debug("Saving the image remotely ...")

        # Debugging
        log.debug("Remote temporary path of image: " + remote_image_path)

        # Save the image remotely
        self.remote.send_python_line(self.label + ".save('" + remote_image_path + "')", show_output=True)

        # Debugging
        log.debug("Downloading the image ...")

        # Download
        self.remote.download(remote_image_path, local_directory)

        # Remove the remote file
        self.remote.remove_file(remote_image_path)

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
        image = self.local_class().from_file(local_path)

        # Remove the local temporary file
        fs.remove_file(local_path)

        # Return the image
        return image

    # -----------------------------------------------------------------

    def convolve(self, kernel, allow_huge=False, fft=True):

        """
        This function ...
        :param kernel:
        :param allow_huge:
        :param fft:
        :return:
        """

        # Inform the user
        log.info("Convolving the remote image ...")

        # SAVE KERNEL LOCALLY

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "kernel.fits")

        # Save the kernel locally
        kernel.save(local_path)

        # UPLOAD KERNEL

        # Remote temp path
        remote_temp_path = self.remote.session_temp_directory

        # Upload the kernel file
        remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
        log.info("Uploading the kernel to " + remote_kernel_path + " ...")
        self.remote.upload(local_path, remote_temp_path, compress=True, show_output=True)

        # Open the kernel remotely
        self.remote.send_python_line("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')", show_output=True)

        # Convolve the image remotely
        self.remote.send_python_line(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ", fft=" + str(fft) + ")", show_output=True, timeout=None)

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

        # Rebin remotely
        self.remote.send_python_line(self.label + ".rebin(reference_wcs, exact=" + str(exact) + ", parallel=" + str(parallel) + ")", timeout=None)

# -----------------------------------------------------------------

class RemoteDataCube(RemoteImage):

    """
    This function ...
    """

    @classmethod
    def from_file(cls, path, wavelength_grid, host_id):

        """
        This function ...
        :param path:
        :param wavelength_grid:
        :param host_id:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Get the remote instance
        remote = prepare_remote(host_id)

        # Remote temp path
        remote_temp_path = remote.session_temp_directory

        # Get file name
        filename = fs.name(path)

        ### UPLOAD DATACUBE

        # Upload the datacube file
        remote_datacube_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the datacube to '" + remote_datacube_path + "' ...")
        remote.upload(path, remote_temp_path, compress=True, show_output=True)

        ### SAVE AND UPLOAD WAVELENGTH GRID

        local_temp_path = tempfile.gettempdir()
        local_wavelength_grid_path = fs.join(local_temp_path, "wavelength_grid.dat")
        log.debug("Saving the wavelength grid locally to '" + local_wavelength_grid_path + "' ...")

        # Save
        wavelength_grid.save(local_wavelength_grid_path)

        # Upload
        remote_wavelength_grid_path = fs.join(remote_temp_path, "wavelength_grid.dat")
        log.info("Uploading the wavelength grid to '" + remote_wavelength_grid_path + "' ...")
        remote.upload(local_wavelength_grid_path, remote_temp_path, show_output=True)

        ###

        # Find label
        label = get_new_label(cls.local_classname(), remote)

        # Create RemoteDataCube instance
        remotedatacube = cls(label, remote)

        # Open the wavelength grid remotely
        log.info("Loading the wavelength grid on the remote host ...")

        # Import the WavelengthGrid class remotely
        remote.import_python_package("WavelengthGrid", from_name="pts.core.simulation.wavelengthgrid")
        remote.send_python_line("wavelength_grid = WavelengthGrid.from_file('" + remote_wavelength_grid_path + "')")

        # Open the frame remotely
        log.info("Loading the datacube on the remote host ...")

        # Actually create the frame remotely
        remote.send_python_line(label + " = " + cls.local_classname() + ".from_file('" + remote_datacube_path + "', wavelength_grid)")

        # Return the remotedatacube instance
        return remotedatacube

    # -----------------------------------------------------------------

    @classmethod
    def from_local(cls, datacube, host_id):

        """
        This function ...
        :param datacube:
        :param host_id:
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "newdatacube_fromlocal.fits")

        # Save the image locally
        datacube.save(local_path)

        # Create the remotedatacube from the locally saved FITS file
        remotedatacube = cls.from_file(local_path, datacube.wavelength_grid, host_id)

        # Remove the local file
        fs.remove_file(local_path)

        # Return the new remotedatacube
        return remotedatacube

    # -----------------------------------------------------------------

    def convolve_with_filters(self, filters):

        """
        This function ...
        :param filters:
        :return:
        """

        # Initialize filter list remotely
        self.remote.send_python_line("filters = []")

        # Reconstruct the list of filters remotely
        for fltr in filters: self.remote.send_python_line("filters.append(Filter.from_string('" + str(fltr) + "'))")

        # Initialize a list with remoteframes
        remoteframes = []

        # Do the convolution remotely
        self.remote.send_python_line("filterconvolvedframes = " + self.label + ".convolve_with_filters(filters)")

        # Create a remoteframe pointing to each of the frames in 'filterconvolvedframes'
        for i in range(len(filters)):

            # Assign a remote label to this result frame
            label_i = get_new_label("Frame", self.remote)

            # Do the assignment remotely
            self.remote.send_python_line(label_i + " = filterconvolvedframes[" + str(i) + "]")

            # Create remoteframe and add it to the list
            remoteframe = RemoteFrame(label_i, self.remote)
            remoteframes.append(remoteframe)

        # Return the list of remoteframes
        return remoteframes

# -----------------------------------------------------------------
