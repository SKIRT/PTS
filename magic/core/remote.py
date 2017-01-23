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

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from .frame import Frame # IMPORTANT THAT THESE ARE IMPORTED !!
from .image import Image # IMPORTANT THAT THESE ARE IMPORTED !!
from .datacube import DataCube # IMPORTANT THAT THESE ARE IMPORTED !!
from ...core.basics.filter import Filter
from ...core.tools import parsing
from ..basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

prepared_sessions = [] # the list of session IDs that have been prepared

# -----------------------------------------------------------------

def prepare_session(session):

    """
    This function ...
    :param session:
    :return:
    """

    # If this host has already been prepared, just return the prepared remote
    if session.session_id in prepared_sessions: return

    # Import
    import_necessary_modules(session)

    # Set log level to debug
    if log.is_debug(): set_debug_log_level(session)

    # Add the host ID to the 'prepared' list
    prepared_sessions.append(session.session_id)

# -----------------------------------------------------------------

def import_necessary_modules(session):

    """
    This function ...
    :param session:
    :return:
    """

    # Inform the user
    log.info("Importing necessary modules ...")

    # Import standard modules
    #session.import_package("tempfile")  ## doesn't work: we seem to have no permissions in this directory on nancy
    session.import_package("urllib")
    session.import_package("numpy", as_name="np")

    # Import the necessary PTS classes and modules
    session.import_package("Frame", from_name="pts.magic.core.frame")
    session.import_package("Image", from_name="pts.magic.core.image")
    session.import_package("DataCube", from_name="pts.magic.core.datacube")
    session.import_package("ConvolutionKernel", from_name="pts.magic.core.kernel")
    session.import_package("CoordinateSystem", from_name="pts.magic.basics.coordinatesystem")
    session.import_package("archive", from_name="pts.core.tools")
    session.import_package("parsing", from_name="pts.core.tools")
    session.import_package("Filter", from_name="pts.core.basics.filter")

# -----------------------------------------------------------------

def set_debug_log_level(session):

    """
    This function ...
    :param session:
    :return:
    """

    # Inform the user
    log.debug("Setting debug logging level remotely ...")

    # Import logging module and setup logger to DEBUG level
    session.import_package("setup_log", from_name="pts.core.tools.logging")
    session.send_python_line("setup_log(level='DEBUG')")

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

def get_labels(classname, session):

    """
    This function ...
    :param classname:
    :param session:
    :return:
    """

    variables = session.variables()

    labels = []

    for variable in variables:
        if variable.startswith(classname.lower()):
            labels.append(variable)

    return labels

# -----------------------------------------------------------------

def get_indices(classname, session):

    """
    This function ...
    :param classname:
    :param session:
    :return:
    """

    return sorted([int(label.split(classname.lower())[1]) for label in get_labels(classname, session)])

# -----------------------------------------------------------------

def get_new_label(classname, session):

    """
    This function ...
    :param classname:
    :param session:
    :return:
    """

    current_indices = get_indices(classname, session)
    return classname.lower() + str(get_first_missing_integer(current_indices))

# -----------------------------------------------------------------

class RemoteFrame(object):

    """
    This class ...
    """

    # Set the default extension
    default_extension = "fits"

    # -----------------------------------------------------------------

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

    def __init__(self, label, session):

        """
        This function ...
        :param label
        :param session:
        """

        self.label = label
        self.session = session

        # The path
        self.path = None

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "__str__()")

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "__repr__()")

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Delete the Frame on the remote from the global namespace
        self.session.remove_variable(self.label)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Check whether floating point value
        if not isinstance(value, float): raise ValueError("Value must be float (is " + str(type(value)) + ")")

        # Multiply remotely
        self.session.send_line(self.label + " *= " + repr(value))

        # Return self
        return self

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Check whether floating point value
        if not isinstance(value, float): raise ValueError("Value must be float (is " + str(type(value)) + ")")

        # Divide remotely
        self.session.send_line(self.label + " /= " + repr(value))

        # Return self
        return self

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__idiv__(value)

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        return Unit(self.session.get_string("str(" + self.label + ".unit)"))

    # -----------------------------------------------------------------

    @unit.setter
    def unit(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        self.session.send_line(self.label + ".unit = '" + str(unit) + "'")

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):

        """
        This function ...
        :return:
        """

        return self.session.evaluate_boolean_expression(self.label + ".has_wcs")

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem(self.session.get_string(self.label + ".wcs.to_header_string()"))

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Set the WCS remotely
        self.session.send_line(self.label + '.wcs = CoordinateSystem("' + wcs.to_header_string() + '")')

    # -----------------------------------------------------------------

    @property
    def filter(self):

        """
        This function ...
        :return:
        """

        # Get the filter
        fltr = Filter(self.session.get_string("str(" + self.label + ".filter)"))

        # Return the filter
        return fltr

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = parsing.quantity(self.session.get_simple_property(self.label, "fwhm"))
        return fwhm

    # -----------------------------------------------------------------

    @fwhm.setter
    def fwhm(self, fwhm):

        """
        This function ...
        :return:
        """

        self.session.send_line(self.label + ".fwhm = parsing.quantity(" + str(fwhm) + ")")

    # -----------------------------------------------------------------

    @classmethod
    def from_url(cls, url, session, index=None, name=None, description=None, plane=None, hdulist_index=None,
                 no_filter=False, fwhm=None, add_meta=True):

        """
        This function ...
        :param url:
        :param session:
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
        prepare_session(session)

        # Inform the user
        log.info("Downloading file " + url + " ...")

        # Local path
        remote_temp_path = session.session_temp_directory
        session.send_line("filename = fs.name(url)")
        session.send_line("local_path = fs.join('" + remote_temp_path + "', filename)")

        # Download
        session.send_line("urllib.urlretrieve(url, local_path)")

        # Get local path (on remote)
        local_path = session.get_string("local_path")

        if local_path.endswith(".fits"): session.send_line("fits_path = local_path")
        else:

            # Inform the user
            log.info("Decompressing kernel file ...")

            # Fits path
            session.send_line("fits_path = fs.join(temp_path, fs.strip_extension(filename))")

            # Decompress the kernel FITS file
            session.send_line("archive.decompress_file(local_path, fits_path)")

            # Remove the compressed file
            session.send_line("fs.remove_file(local_path)")

        # Find label
        label = get_new_label(cls.local_classname(), session)

        # Create RemoteFrame instance
        remoteframe = cls(label, session)

        # Actually create the frame remotely
        session.send_line(label + " = " + cls.local_classname() + ".from_file(fits_path)")

        # Return the remoteframe instance
        return remoteframe

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, session, index=None, name=None, description=None, plane=None, hdulist_index=None,
                  no_filter=False, fwhm=None, add_meta=True):

        """
        This function ...
        :param path:
        :param session:
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

        # Prepare session if necessary
        prepare_session(session)

        # Remote temp path
        remote_temp_path = session.session_temp_directory

        # Get file name
        filename = fs.name(path)

        # Upload the frame file
        remote_frame_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the frame to " + remote_frame_path + " ...")
        session.remote.upload(path, remote_temp_path, compress=True, show_output=True)

        # Find label
        label = get_new_label(cls.local_classname(), session)

        # Create RemoteFrame instance
        remoteframe = cls(label, session)

        # Open the frame remotely
        log.info("Loading the frame on the remote host ...")

        # Actually create the frame remotely
        session.send_line(label + " = " + cls.local_classname() + ".from_file('" + remote_frame_path + "')")

        # Set the path
        remoteframe.path = path

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
        label = get_new_label(self.local_classname(), self.session)

        # Copy the frame remotely
        self.session.send_line(label + " = " + self.label + ".copy()")

        # Create new RemoteFrame instance
        newremoteframe = self.__class__(label, self.session)

        # Return the new remoteframe instance
        return newremoteframe

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "xsize")

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "ysize")

    # -----------------------------------------------------------------

    @property
    def fwhm_pix(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "fwhm_pix")

    # -----------------------------------------------------------------

    def sum(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "sum()")

    # -----------------------------------------------------------------

    def normalize(self, to=1.0):

        """
        This function ...
        :param to:
        :return:
        """

        # Call the function on the remote
        self.session.send_line(self.label + ".normalize(" + str(to) + ")")

    # -----------------------------------------------------------------

    def is_constant(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "is_constant()")

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
        kernel.saveto(local_path)

        # UPLOAD KERNEL

        # Remote temp path
        remote_temp_path = self.session.session_temp_directory

        # Upload the kernel file
        remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
        log.info("Uploading the kernel to " + remote_kernel_path + " ...")
        self.session.remote.upload(local_path, remote_temp_path, compress=True, show_output=True)

        # Open the kernel remotely
        self.session.send_line("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')")

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
        self.session.send_line(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ", fft=" + str(fft) + ")", show_output=True, timeout=None)

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

    def rebin(self, reference_wcs, exact=False, parallel=True):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :return:
        """

        # Check whether the remote frame has a WCS
        if not self.has_wcs: raise RuntimeError("Cannot rebin a frame without coordinate system")

        # Upload the WCS
        self.session.send_line('reference_wcs = CoordinateSystem("' + reference_wcs.to_header_string() + '")')

        # Create remote frame label for the footprint
        footprint_label = get_new_label("Frame", self.session)

        # Rebin, get the footprint
        self.session.send_line(footprint_label + " = " + self.label + ".rebin(reference_wcs, exact=" + str(exact) + ", parallel=" + str(parallel) + ")", timeout=None, show_output=True)

        # Create a new remoteframe instance for the footprint
        remote_footprint_frame = RemoteFrame(footprint_label, self.session)

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
        self.session.send_line(self.label + ".crop(" + str(x_min) + ", " + str(x_max) + ", " + str(y_min) + ", " + str(y_max))

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
        self.session.send_line(self.label + ".pad(" + str(nx) + ", " + str(ny) + ")")

    # -----------------------------------------------------------------

    def unpad(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        self.session.send_line(self.label + ".unpad(" + str(nx) + ", " + str(ny) + ")")

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

        self.session.send_line(self.label + ".downsample(" + str(factor) + ", " + str(order) + ")")

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

        self.session.send_line(self.label + ".upsample(" + str(factor) + ", " + str(integers) + ")")

    # -----------------------------------------------------------------

    def fill(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line(self.label + ".fill(" + str(value) + ")")

    # -----------------------------------------------------------------

    def replace_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line(self.label + ".replace_nans(" + repr(value) + ")")

    # -----------------------------------------------------------------

    def replace_infs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line(self.label + ".replace_infs(" + repr(value) + ")")

    # -----------------------------------------------------------------

    def replace_negatives(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line(self.label + ".replace_negatives(" + repr(value) + ")")

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
        frame.saveto(local_path)

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
        self.saveto(local_path)

        # Create the local frame
        frame = self.local_class().from_file(local_path)

        # Remove the local temporary file
        fs.remove_file(local_path)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the remote frame ...")

        # Check whether the path is valid
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

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
        remote_temp_path = self.session.session_temp_directory

        # Remote file path
        remote_file_path = fs.join(remote_temp_path, filename)

        # Debugging
        log.debug("Saving the frame remotely ...")

        # Save the frame remotely
        self.session.send_line(self.label + ".saveto('" + remote_file_path + "')")

        # Debugging
        log.debug("Downloading the frame ...")

        # Download
        self.session.remote.download(remote_file_path, local_directory, compress=True, show_output=True)

        # Remove the remote file
        self.session.remove_file(remote_file_path)

        # Update the path
        self.path = path

# -----------------------------------------------------------------

class RemoteImage(object):

    """
    This class ...
    """

    # Set the default extension
    default_extension = "fits"

    # -----------------------------------------------------------------

    @classmethod
    def local_class(cls):

        """
        This function ...
        :return:
        """

        classname = cls.local_classname()
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

    def __init__(self, label, session):

        """
        The constructor ...
        :param label:
        :param session:
        """

        self.label = label
        self.session = session

        # The path
        self.path = None

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "__str__()")

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "__repr__()")

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Delete the Image on the remote from the global namespace
        self.session.remove_variable(self.label)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Check whether floating point value
        if not isinstance(value, float): raise ValueError("Value must be float (is " + str(type(value)) + ")")

        # Multiply remotely
        self.session.send_line(self.label + " *= " + repr(value))

        # Return self
        return self

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Check whether floating point value
        if not isinstance(value, float): raise ValueError("Value must be float (is " + str(type(value)) + ")")

        # Divide remotely
        self.session.send_line(self.label + " /= " + repr(value))

        # Return self
        return self

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__idiv__(value)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, session):

        """
        This function ...
        :param path:
        :param session:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Get the remote instance
        prepare_session(session)

        # Remote temp path
        remote_temp_path = session.session_temp_directory

        # Get file name
        filename = fs.name(path)

        # Upload the image file
        remote_image_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the image to '" + remote_image_path + "' ...")
        session.remote.upload(path, remote_temp_path, compress=True, show_output=True)

        # Find label
        label = get_new_label(cls.local_classname(), session)

        # Create RemoteImage instance
        remoteimage = cls(label, session)

        # Open the frame remotely
        log.info("Loading the image on the remote host ...")

        # Actually create the frame remotely
        session.send_line(label + " = " + cls.local_classname() + ".from_file('" + remote_image_path + "')")

        # Set the path
        remoteimage.path = path

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
        label = get_new_label(self.local_classname(), self.session)

        # Copy the frame remotely
        self.session.send_line(label + " = " + self.label + ".copy()")

        # Create new RemoteImage instance
        newremoteimage = self.__class__(label, self.session)

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
        label = get_new_label("Frame", self.session)

        # Reference the primary frame to the new variable
        self.session.send_python_line(label + " = " + self.label + ".primary")

        # Create a new RemoteFrame instance
        newremoteframe = RemoteFrame(label, self.session) # We assume the frame of the remote image is not of a type derived from Frame

        # Return the new remoteframe instance
        return newremoteframe

    # -----------------------------------------------------------------

    @property
    def has_frames(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simpl_property(self.label, "has_frames")

    # -----------------------------------------------------------------

    @property
    def nframes(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "nframes")

    # -----------------------------------------------------------------

    @property
    def nmasks(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "nmasks")

    # -----------------------------------------------------------------

    @property
    def nregions(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "nregions")

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "shape")

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        return Unit(self.session.get_string("str(" + self.label + ".unit)"))

    # -----------------------------------------------------------------

    @unit.setter
    def unit(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        self.session.send_python_line(self.label + ".unit = '" + str(unit) + "'")

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):

        """
        This function ...
        :return:
        """

        return self.session.evaluate_boolean_expression(self.label + ".has_wcs")

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem(self.session.get_string(self.label + ".wcs.to_header_string()"))

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Set the WCS remotely
        self.session.send_line(self.label + '.wcs = CoordinateSystem("' + wcs.to_header_string() + '")')

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = parsing.quantity(self.session.get_simple_property(self.label, "fwhm"))
        return fwhm

    # -----------------------------------------------------------------

    @fwhm.setter
    def fwhm(self, fwhm):

        """
        This function ...
        :return:
        """

        self.session.send_line(self.label + ".fwhm = parsing.quantity(" + str(fwhm) + ")")

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the remote image ...")

        # Check whether the path is valid
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

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
        remote_temp_path = self.session.session_temp_directory

        # Determine remote path for the image
        remote_image_path = fs.join(remote_temp_path, filename)

        # Debugging
        log.debug("Saving the image remotely ...")

        # Debugging
        log.debug("Remote temporary path of image: " + remote_image_path)

        # Save the image remotely
        self.session.send_line(self.label + ".saveto('" + remote_image_path + "')", show_output=True)

        # Debugging
        log.debug("Downloading the image ...")

        # Download
        self.session.remote.download(remote_image_path, local_directory, compress=True, show_output=True)

        # Remove the remote file
        self.session.remove_file(remote_image_path)

        # Update the path
        self.path = path

    # -----------------------------------------------------------------

    @classmethod
    def from_local(cls, image, session):

        """
        This function ...
        :param image:
        :param session:
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "newimage_fromlocal.fits")

        # Save the image locally
        image.saveto(local_path)

        # Create the remoteimage from the locally saved image
        remoteimage = cls.from_file(local_path, session)

        # Remove the local file
        fs.remove_file(local_path)

        # Set the path
        remoteimage.path = image.path

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
        self.saveto(local_path)

        # Create the image
        image = self.local_class().from_file(local_path)

        # Remove the local temporary file
        fs.remove_file(local_path)

        # Set the path
        image.path = self.path

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
        kernel.saveto(local_path)

        # UPLOAD KERNEL

        # Remote temp path
        remote_temp_path = self.session.session_temp_directory

        # Upload the kernel file
        remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
        log.info("Uploading the kernel to " + remote_kernel_path + " ...")
        self.session.remote.upload(local_path, remote_temp_path, compress=True, show_output=True)

        # Open the kernel remotely
        self.session.send_line("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')", show_output=True)

        # Convolve the image remotely
        self.session.send_line(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ", fft=" + str(fft) + ")", show_output=True, timeout=None)

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=False, parallel=True):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :return:
        """

        # Check whether the remote image has a WCS
        if not self.has_wcs: raise RuntimeError("Cannot rebin a frame without coordinate system")

        # Upload the WCS
        self.session.send_line('reference_wcs = CoordinateSystem("' + reference_wcs.to_header_string() + '")')

        # Rebin remotely
        self.session.send_line(self.label + ".rebin(reference_wcs, exact=" + str(exact) + ", parallel=" + str(parallel) + ")", show_output=True, timeout=None)

# -----------------------------------------------------------------

class RemoteDataCube(RemoteImage):

    """
    This function ...
    """

    # Set the default extension
    default_extension = "fits"

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, wavelength_grid, session):

        """
        This function ...
        :param path:
        :param wavelength_grid:
        :param session:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Get the remote instance
        prepare_session(session)

        # Remote temp path
        remote_temp_path = session.session_temp_directory

        # Get file name
        filename = fs.name(path)

        ### UPLOAD DATACUBE

        # Upload the datacube file
        remote_datacube_path = fs.join(remote_temp_path, filename)
        log.info("Uploading the datacube to '" + remote_datacube_path + "' ...")
        session.remote.upload(path, remote_temp_path, compress=True, show_output=True)

        ### SAVE AND UPLOAD WAVELENGTH GRID

        local_temp_path = tempfile.gettempdir()
        local_wavelength_grid_path = fs.join(local_temp_path, "wavelength_grid.dat")
        log.debug("Saving the wavelength grid locally to '" + local_wavelength_grid_path + "' ...")

        # Save
        wavelength_grid.saveto(local_wavelength_grid_path)

        # Upload
        remote_wavelength_grid_path = fs.join(remote_temp_path, "wavelength_grid.dat")
        log.info("Uploading the wavelength grid to '" + remote_wavelength_grid_path + "' ...")
        session.remote.upload(local_wavelength_grid_path, remote_temp_path, show_output=True)

        ###

        # Find label
        label = get_new_label(cls.local_classname(), session)

        # Create RemoteDataCube instance
        remotedatacube = cls(label, session)

        # Open the wavelength grid remotely
        log.info("Loading the wavelength grid on the remote host ...")

        # Import the WavelengthGrid class remotely
        session.import_package("WavelengthGrid", from_name="pts.core.simulation.wavelengthgrid")
        session.send_line("wavelength_grid = WavelengthGrid.from_file('" + remote_wavelength_grid_path + "')")

        # Open the frame remotely
        log.info("Loading the datacube on the remote host ...")

        # Actually create the frame remotely
        session.send_line(label + " = " + cls.local_classname() + ".from_file('" + remote_datacube_path + "', wavelength_grid)", show_output=True)

        # Return the remotedatacube instance
        return remotedatacube

    # -----------------------------------------------------------------

    @classmethod
    def from_local(cls, datacube, session):

        """
        This function ...
        :param datacube:
        :param session:
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "newdatacube_fromlocal.fits")

        # Save the image locally
        datacube.saveto(local_path)

        # Create the remotedatacube from the locally saved FITS file
        remotedatacube = cls.from_file(local_path, datacube.wavelength_grid, session)

        # Remove the local file
        fs.remove_file(local_path)

        # Return the new remotedatacube
        return remotedatacube

    # -----------------------------------------------------------------

    def convolve_with_filters(self, filters, nprocesses=8):

        """
        This function ...
        :param filters:
        :param nprocesses:
        :return:
        """

        # Initialize filter list remotely
        self.session.send_line("filters = []")

        # Reconstruct the list of filters remotely
        for fltr in filters: self.session.send_line("filters.append(Filter('" + str(fltr) + "'))")

        # Initialize a list with remoteframes
        remoteframes = []

        # Do the convolution remotely
        self.session.send_line("filterconvolvedframes = " + self.label + ".convolve_with_filters(filters, nprocesses=" + str(nprocesses) + ")", timeout=None, show_output=True)

        # Create a remoteframe pointing to each of the frames in 'filterconvolvedframes'
        for i in range(len(filters)):

            # Assign a remote label to this result frame
            label_i = get_new_label("Frame", self.session)

            # Do the assignment remotely
            self.session.send_line(label_i + " = filterconvolvedframes[" + str(i) + "]")

            # Create remoteframe and add it to the list
            remoteframe = RemoteFrame(label_i, self.session)
            remoteframes.append(remoteframe)

        # Return the list of remoteframes
        return remoteframes

    # -----------------------------------------------------------------

    def to_wavelength_density(self, new_unit, wavelength_unit):

        """
        This function ...
        :param new_unit:
        :param wavelength_unit:
        :return:
        """

        # Convert to wavelength density remotely
        self.session.send_line(self.label + ".to_wavelength_density('" + new_unit + "', '" + wavelength_unit + "')")

# -----------------------------------------------------------------
