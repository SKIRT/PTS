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
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from .frame import Frame # IMPORTANT THAT THESE ARE IMPORTED !!
from .image import Image # IMPORTANT THAT THESE ARE IMPORTED !!
from .datacube import DataCube # IMPORTANT THAT THESE ARE IMPORTED !!
from .datacube import needs_spectral_convolution
from ...core.filter.filter import parse_filter
from ...core.tools import parsing
from ..basics.coordinatesystem import CoordinateSystem
from ...core.units.parsing import parse_unit as u
from ...core.tools import types
from ..basics.pixelscale import Pixelscale, PhysicalPixelscale

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

    # Import standard modules
    session.import_package_update("urllib", show_output=log.is_debug())
    session.import_package_update("numpy", as_name="np", show_output=log.is_debug())

    # Import the necessary PTS classes and modules
    session.import_package_update("Frame", from_name="pts.magic.core.frame", show_output=log.is_debug())
    session.import_package_update("Image", from_name="pts.magic.core.image", show_output=log.is_debug())
    session.import_package_update("DataCube", from_name="pts.magic.core.datacube", show_output=log.is_debug())
    session.import_package_update("ConvolutionKernel", from_name="pts.magic.core.kernel", show_output=log.is_debug())
    session.import_package_update("CoordinateSystem", from_name="pts.magic.basics.coordinatesystem", show_output=log.is_debug())
    session.import_package_update("archive", from_name="pts.core.tools", show_output=log.is_debug())
    session.import_package_update("parsing", from_name="pts.core.tools", show_output=log.is_debug())
    session.import_package_update("parse_filter", from_name="pts.core.filter.filter", show_output=log.is_debug())
    session.import_package_update("BroadBandFilter", from_name="pts.core.filter.broad", show_output=log.is_debug())
    session.import_package_update("tostr", from_name="pts.core.tools.stringify", show_output=log.is_debug())

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
    session.import_package("setup_log", from_name="pts.core.basics.log")
    session.send_line_and_raise("setup_log(level='DEBUG')")

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

def get_filter(frame_path, session):

    """
    Ths function allows getting the filter of a frame without loading the entire frame
    :param frame_path:
    :param session:
    :return:
    """

    name = get_filter_name(frame_path, session)
    if name is None: return None
    else: return parse_filter(name)

# -----------------------------------------------------------------

def get_filter_name(frame_path, session):

    """
    This function ...
    :param frame_path:
    :param session:
    :return:
    """

    # header = fits.getheader(frame_path)
    session.import_package("getheader", from_name="astropy.io.fits")
    session.send_line_and_raise("header = getheader('" + frame_path + "')")
    session.import_package("get_filter", from_name="pts.magic.tools.headers")
    name = fs.name(frame_path[:-5])
    # fltr = headers.get_filter(fs.name(frame_path[:-5]), header)
    return session.get_simple_variable("str(get_filter('" + name + "', header))")

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
        self.session.send_line_and_raise(self.label + " *= " + repr(value))

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
        self.session.send_line_and_raise(self.label + " /= " + repr(value))

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

        return u(self.session.get_string("str(" + self.label + ".unit)"))

    # -----------------------------------------------------------------

    @unit.setter
    def unit(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        self.session.send_line_and_raise(self.label + ".unit = '" + str(unit) + "'")

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

        return CoordinateSystem.from_header_string(self.session.get_string(self.label + ".wcs.to_header_string()"))

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Set the WCS remotely
        self.session.send_line_and_raise(self.label + '.wcs = CoordinateSystem.from_header_string("' + wcs.to_header_string() + '")')

    # -----------------------------------------------------------------

    @property
    def filter(self):

        """
        This function ...
        :return:
        """

        # Get the filter
        fltr = parse_filter(self.session.get_string("str(" + self.label + ".filter)"))

        # Return the filter
        return fltr

    # -----------------------------------------------------------------

    @filter.setter
    def filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Set the filter
        self.session.send_line_and_raise(self.label + ".filter = '" + str(fltr) + "'")

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

        self.session.send_line_and_raise(self.label + ".fwhm = parsing.quantity(" + str(fwhm) + ")")

    # -----------------------------------------------------------------

    @property
    def has_pixelscale(self):

        """
        Thisfunction ...
        :return:
        """

        string = self.session.get_string("tostr(" + self.label + ".pixelscale.x)")
        return string != "None"

    # -----------------------------------------------------------------

    @property
    def x_pixelscale(self):

        """
        This function ...
        :return:
        """

        string = self.session.get_string("tostr(" + self.label + ".pixelscale.x)")
        return parsing.angle_or_length_quantity(string)

    # -----------------------------------------------------------------

    @property
    def y_pixelscale(self):

        """
        This function ...
        :return:
        """

        string = self.session.get_string("tostr(" + self.label + ".pixelscale.y)")
        return parsing.angle_or_length_quantity(string)

    # -----------------------------------------------------------------

    @property
    def has_angular_pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: raise ValueError("No pixelscale")
        return types.is_angle(self.x_pixelscale)

    # -----------------------------------------------------------------

    @property
    def has_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: raise ValueError("No pixelscale")
        return types.is_length_quantity(self.x_pixelscale)

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: return None
        else:
            if self.has_angular_pixelscale: return Pixelscale(x=self.x_pixelscale, y=self.y_pixelscale)
            elif self.has_physical_pixelscale: return PhysicalPixelscale(x=self.x_pixelscale, y=self.y_pixelscale)
            else: raise ValueError("Uknown pixelscale type")

    # -----------------------------------------------------------------

    @property
    def average_pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: return None
        else:
            string = self.session.get_string("tostr(" + self.label + ".average_pixelscale)")
            return parsing.angle_or_length_quantity(string)

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
        session.send_line_and_raise("filename = fs.name(url)")
        session.send_line_and_raise("local_path = fs.join('" + remote_temp_path + "', filename)")

        # Download
        session.send_line_and_raise("urllib.urlretrieve(url, local_path)")

        # Get local path (on remote)
        local_path = session.get_string("local_path")

        if local_path.endswith(".fits"): session.send_line_and_raise("fits_path = local_path")
        else:

            # Inform the user
            log.info("Decompressing kernel file ...")

            # Fits path
            session.send_line_and_raise("fits_path = fs.join(temp_path, fs.strip_extension(filename))")

            # Decompress the kernel FITS file
            session.send_line_and_raise("archive.decompress_file(local_path, fits_path)")

            # Remove the compressed file
            session.send_line_and_raise("fs.remove_file(local_path)")

        # Find label
        label = get_new_label(cls.local_classname(), session)

        # Create RemoteFrame instance
        remoteframe = cls(label, session)

        # Actually create the frame remotely
        session.send_line_and_raise(label + " = " + cls.local_classname() + ".from_file(fits_path)")

        # Return the remoteframe instance
        return remoteframe

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, session):

        """
        This function ...
        :param path:
        :param session:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file '" + path + "' ...")

        # Prepare session if necessary
        prepare_session(session)

        # Get file name
        #filename = fs.name(path)

        # Find label
        label = get_new_label(cls.local_classname(), session)

        # Create RemoteFrame instance
        remoteframe = cls(label, session)

        # Open the frame remotely
        log.info("Loading the frame on the remote host ...")

        # Actually create the frame remotely
        session.send_line_and_raise(label + " = " + cls.local_classname() + ".from_file('" + path + "')")

        # Set the path
        remoteframe.path = path

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
        log.info("Reading in file '" + path + "' ...")

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
        session.send_line_and_raise(label + " = " + cls.local_classname() + ".from_file('" + remote_frame_path + "')")

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
        self.session.send_line_and_raise(label + " = " + self.label + ".copy()")

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

    @property
    def sigma(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "sigma")

    # -----------------------------------------------------------------

    @property
    def sigma_pix(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "sigma_pix")

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
        :return:
        """

        from ...core.units.unit import PhotometricUnit
        from ...core.tools import types
        from ...core.tools.stringify import tostr

        # Determine distance string
        if distance is not None:
            distance_string = tostr(distance)
            self.session.import_package("parse_quantity", from_name="pts.core.units.parsing")
            parse_distance_string = "parse_quantity('" + distance_string + "')"
        else: parse_distance_string = "None"

        # Photometric unit
        if isinstance(to_unit, PhotometricUnit):

            # Get unit string
            to_unit_string = tostr(to_unit)

            # Density unit
            if to_unit.density:
                if not density and density_strict: raise ValueError("Density = False and density_strict = True is incompatible with the passed unit")
                density = True
                density_strict = True

            # Not density unit
            else:
                if density and density_strict: raise ValueError("Density = True and density_strict = True is incompatible with the passed unit")
                density = False
                density_strict = True

            # Brightness unit
            if to_unit.brightness:
                if not brightness and brightness_strict: raise ValueError("Brightness = False and brightness_strict = True is incompatible with the passed unit")
                brightness = True
                brightness_strict = True

            # Not a brightness unit
            else:
                if brightness and brightness_strict: raise ValueError("Brightness = True and brightness_strict = True is incompatible with the passed unit")
                brightness = False
                brightness_strict = True

        # String
        elif types.is_string_type(to_unit): to_unit_string = to_unit

        # Invalid
        else: raise ValueError("Invalid value for 'to_unit': must be photometric unit or string")

        # Send command
        self.session.send_line_and_raise(self.label + ".convert_to('" + to_unit_string + "', distance=" + parse_distance_string + ", density=" + tostr(density) + ", brightness=" + tostr(brightness) + ", density_strict=" + tostr(density_strict) + ", brightness_strict=" + tostr(brightness_strict) + ")")

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
        self.session.send_line_and_raise(self.label + ".normalize(" + str(to) + ")")

    # -----------------------------------------------------------------

    @property
    def is_constant(self):

        """
        This function ...
        :return:
        """

        return self.session.get_simple_property(self.label, "is_constant")

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

    def convolve(self, kernel, allow_huge=True, fft=True):

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
        if self.is_constant:
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
        self.session.send_line_and_raise("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')")

        # Prepare the kernel if necessary
        #self.remote.send_line("if not kernel.prepared: kernel.prepare_for(" + self.label + ")")

        # Check where the NaNs are at
        #self.remote.send_line("nans_mask = np.isnan(" + self.label + "._data)")

        # Assert that the kernel is normalized
        #self.remote.send_line("assert kernel.normalized")

        # Do the convolution
        #self.remote. .... ETC


        # OR JUST: (does constant check again and creates copy for the prepared kernel, but... much more maintainable in this way)

        # Do the convolution on the remote frame
        self.session.send_line_and_raise(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ", fft=" + str(fft) + ")", show_output=True, timeout=None)

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
        self.session.send_line_and_raise('reference_wcs = CoordinateSystem.from_header_string("' + reference_wcs.to_header_string() + '")')

        # Create remote frame label for the footprint
        footprint_label = get_new_label("Frame", self.session)

        # Rebin, get the footprint
        self.session.send_line_and_raise(footprint_label + " = " + self.label + ".rebin(reference_wcs, exact=" + str(exact) + ", parallel=" + str(parallel) + ")", timeout=None, show_output=True)

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
        self.session.send_line_and_raise(self.label + ".crop(" + str(x_min) + ", " + str(x_max) + ", " + str(y_min) + ", " + str(y_max))

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
        self.session.send_line_and_raise(self.label + ".pad(" + str(nx) + ", " + str(ny) + ")")

    # -----------------------------------------------------------------

    def unpad(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        self.session.send_line_and_raise(self.label + ".unpad(" + str(nx) + ", " + str(ny) + ")")

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

        self.session.send_line_and_raise(self.label + ".downsample(" + str(factor) + ", " + str(order) + ")")

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

        self.session.send_line_and_raise(self.label + ".upsample(" + str(factor) + ", " + str(integers) + ")")

    # -----------------------------------------------------------------

    def fill(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line_and_raise(self.label + ".fill(" + str(value) + ")")

    # -----------------------------------------------------------------

    def replace_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line_and_raise(self.label + ".replace_nans(" + repr(value) + ")")

    # -----------------------------------------------------------------

    def replace_infs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line_and_raise(self.label + ".replace_infs(" + repr(value) + ")")

    # -----------------------------------------------------------------

    def replace_negatives(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.session.send_line_and_raise(self.label + ".replace_negatives(" + repr(value) + ")")

    # -----------------------------------------------------------------

    @classmethod
    def from_local(cls, frame, session):

        """
        This function ...
        :param frame:
        :param session:
        :return:
        """

        # Determine temporary directory path
        temp_path = tempfile.gettempdir()

        # Determine local FITS path
        local_path = fs.join(temp_path, "newframe_fromlocal.fits")

        # Save the frame locally
        frame.saveto(local_path)

        # Create the remoteframe from the locally saved frame
        remoteframe = cls.from_file(local_path, session)

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

        # SAVE REMOTELY
        self.saveto_remote(remote_file_path)

        # Debugging
        log.debug("Downloading the frame ...")

        # Download
        self.session.remote.download(remote_file_path, local_directory, compress=True, show_output=True)

        # Remove the remote file
        self.session.remove_file(remote_file_path)

        # Update the path
        self.path = path

    # -----------------------------------------------------------------

    def saveto_remote(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Saving the frame remotely ...")

        # Save the frame remotely
        self.session.send_line_and_raise(self.label + ".saveto('" + path + "')", show_output=True)

    # -----------------------------------------------------------------

    def saveto_png(self, path, interval="pts", scale="log", alpha="absolute", peak_alpha=1., colours="red"):

        """
        This function ...
        :param path:
        :param interval:
        :param scale:
        :param alpha:
        :param peak_alpha:
        :param colours:
        :return:
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
        self.session.send_line_and_raise(self.label + ".saveto_png('" + remote_file_path + "', interval='" + interval + "', scale='" + scale + "', alpha='" + str(alpha) + "', peak_alpha=" + str(peak_alpha) + ", colours='" + colours + "')", show_output=True)

        # Debugging
        log.debug("Downloading the image ...")

        # Download
        self.session.remote.download(remote_file_path, local_directory, compress=False, show_output=True)

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
        self.session.send_line_and_raise(self.label + " *= " + repr(value))

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
        self.session.send_line_and_raise(self.label + " /= " + repr(value))

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
        log.info("Reading in file '" + path + "' ...")

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
        session.send_line_and_raise(label + " = " + cls.local_classname() + ".from_file('" + remote_image_path + "')")

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
        self.session.send_line_and_raise(label + " = " + self.label + ".copy()")

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
        self.session.send_line_and_raise(label + " = " + self.label + ".primary")

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

        return u(self.session.get_string("str(" + self.label + ".unit)"))

    # -----------------------------------------------------------------

    @unit.setter
    def unit(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        self.session.send_line_and_raise(self.label + ".unit = '" + str(unit) + "'")

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

        return CoordinateSystem.from_header_string(self.session.get_string(self.label + ".wcs.to_header_string()"))

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Set the WCS remotely
        self.session.send_line_and_raise(self.label + '.wcs = CoordinateSystem.from_header_string("' + wcs.to_header_string() + '")')

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

        self.session.send_line_and_raise(self.label + ".fwhm = parsing.quantity(" + str(fwhm) + ")")

    # -----------------------------------------------------------------

    def convert_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False, brightness_strict=False):

        """
        This function ...
        """

        from ...core.units.unit import PhotometricUnit
        from ...core.tools import types
        from ...core.tools.stringify import tostr

        # Determine distance string
        if distance is not None:
            distance_string = tostr(distance)
            self.session.import_package("parse_quantity", from_name="pts.core.units.parsing")
            parse_distance_string = "parse_quantity('" + distance_string + "')"
        else: parse_distance_string = "None"

        # Photometric unit
        if isinstance(to_unit, PhotometricUnit):

            # Get unit string
            to_unit_string = tostr(to_unit)

            # Density unit
            if to_unit.density:
                if not density and density_strict: raise ValueError("Density = False and density_strict = True is incompatible with the passed unit")
                density = True
                density_strict = True

            # Not density unit
            else:
                if density and density_strict: raise ValueError("Density = True and density_strict = True is incompatible with the passed unit")
                density = False
                density_strict = True

            # Brightness unit
            if to_unit.brightness:
                if not brightness and brightness_strict: raise ValueError("Brightness = False and brightness_strict = True is incompatible with the passed unit")
                brightness = True
                brightness_strict = True

            # Not a brightness unit
            else:
                if brightness and brightness_strict: raise ValueError("Brightness = True and brightness_strict = True is incompatible with the passed unit")
                brightness = False
                brightness_strict = True

        # String
        elif types.is_string_type(to_unit): to_unit_string = to_unit

        # Invalid
        else: raise ValueError("Invalid value for 'to_unit': must be photometric unit or string")

        # Send command
        self.session.send_line_and_raise(self.label + ".convert_to('" + to_unit_string + "', distance=" + parse_distance_string + ", density=" + tostr(density) + ", brightness=" + tostr(brightness) + ", density_strict=" + tostr(density_strict) + ", brightness_strict=" + tostr(brightness_strict) + ")")

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
        log.debug("Remote temporary path of image: " + remote_image_path)

        # SAVE REMOTELY
        self.saveto_remote(remote_image_path)

        # Debugging
        log.debug("Downloading the image ...")

        # Download
        self.session.remote.download(remote_image_path, local_directory, compress=True, show_output=True)

        # Remove the remote file
        self.session.remove_file(remote_image_path)

        # Update the path
        self.path = path

    # -----------------------------------------------------------------

    def saveto_remote(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Saving the image remotely ...")

        # Save the image remotely
        self.session.send_line_and_raise(self.label + ".saveto('" + path + "')", show_output=True)

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

    def convolve(self, kernel, allow_huge=True, fft=True):

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
        self.session.send_line_and_raise("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')", show_output=True)

        # Convolve the image remotely
        self.session.send_line_and_raise(self.label + ".convolve(kernel, allow_huge=" + str(allow_huge) + ", fft=" + str(fft) + ")", show_output=True, timeout=None)

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
        self.session.send_line_and_raise('reference_wcs = CoordinateSystem.from_header_string("' + reference_wcs.to_header_string() + '")')

        # Rebin remotely
        self.session.send_line_and_raise(self.label + ".rebin(reference_wcs, exact=" + str(exact) + ", parallel=" + str(parallel) + ")", show_output=True, timeout=None)

# -----------------------------------------------------------------

class RemoteDataCube(RemoteImage):

    """
    This function ...
    """

    # Set the default extension
    default_extension = "fits"

    # -----------------------------------------------------------------

    @staticmethod
    def find_in_previous_session(path, session):

        """
        This function ...
        :param path:
        :param session:
        :return:
        """

        filename = fs.name(path)

        # FIRST CHECK WHETHER THE FILE WAS ALREADY UPLOADED IN A PREVIOUS SESSION
        # Loop over the previous sessions in reversed order
        for session_path in session.previous_temp_directories:

            # Determine the potential path to the datacube
            remote_datacube_path = fs.join(session_path, filename)

            # Determine the potential path to the wavelength grid
            remote_wavelength_grid_path = fs.join(session_path, "wavelength_grid.dat")

            # Check whether the file exists and the same as the local file
            if session.remote.exists_and_equal_to_local(remote_datacube_path, path) and session.remote.is_file(remote_wavelength_grid_path):

                # Inform the user
                log.success("The datacube is found to be already uploaded in previous session '" + fs.name(session_path) + "': not uploading ...")

                # Return the paths
                return remote_datacube_path, remote_wavelength_grid_path

        # Nothing found
        return None, None

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
        log.info("Reading in file '" + path + "' ...")

        # Get file name
        filename = fs.name(path)

        # Get the remote instance
        prepare_session(session)

        # Find in previous session
        remote_datacube_path, remote_wavelength_grid_path = cls.find_in_previous_session(path, session)
        #print(remote_datacube_path, remote_wavelength_grid_path)

        # Not found
        if remote_datacube_path is None:

            # Remote temp path
            remote_temp_path = session.session_temp_directory

            ## CHECK WHETHER THE FILE IS VALID BEFORE UPLOADING!!
            from . import fits
            if not fits.is_valid(path): raise fits.DamagedFITSFileError("Local FITS file is damaged", path=path)

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

        # Find new label
        label = get_new_label(cls.local_classname(), session)

        # Create RemoteDataCube instance
        remotedatacube = cls(label, session)

        # Open the wavelength grid remotely
        log.info("Loading the wavelength grid on the remote host ...")

        # Import the WavelengthGrid class remotely
        session.import_package("WavelengthGrid", from_name="pts.core.simulation.wavelengthgrid")
        session.send_line_and_raise("wavelength_grid = WavelengthGrid.from_file('" + remote_wavelength_grid_path + "')")

        # Open the frame remotely
        log.info("Loading the datacube on the remote host ...")

        # Actually create the frame remotely
        session.send_line_and_raise(label + " = " + cls.local_classname() + ".from_file('" + remote_datacube_path + "', wavelength_grid)", show_output=True)

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

    def get_frame_index_for_wavelength(self, wavelength):

        """
        This function ...
        :return:
        """

        # Import the unit parsing function rmeotely
        self.session.import_package("parse_unit", from_name="pts.core.units.parsing")
        self.session.import_package("parse_quantity", from_name="pts.core.units.parsing")
        #self.session.import_package("tostr", from_name="pts.core.tools.stringify")
        from pts.core.tools.stringify import tostr
        return self.session.get_simple_property(self.label, "get_frame_index_for_wavelength(parse_quantity('" + tostr(wavelength) + "'))")

    # -----------------------------------------------------------------

    def frames_for_filters(self, filters, convolve=False, nprocesses=8, check_previous_sessions=False):

        """
        This function ...
        :param filters:
        :param convolve:
        :param nprocesses:
        :param check_previous_sessions:
        :return:
        """

        # Inform the user
        log.info("Getting frames for " + str(len(filters)) + " different filters ...")

        remoteframes = []
        for_convolution = []

        # Loop over the filters
        for fltr in filters:

            # Needs spectral convolution
            if needs_spectral_convolution(fltr, convolve):

                # Debugging
                log.debug("The frame for the " + str(fltr) + " filter will be calculated by convolving spectrally")

                # Add to list
                for_convolution.append(fltr)

                # Add placeholder
                remoteframes.append(None)

            # No spectral convolution for this filter
            else:

                # Debugging
                log.debug("Getting the frame for the " + str(fltr) + " filter ...")

                # Get the index of the wavelength closest to that of the filter
                index = self.get_frame_index_for_wavelength(fltr.pivot)

                # Assign a remote label to this result frame
                label_i = get_new_label("Frame", self.session)

                # Do the assignment remotely
                self.session.send_line_and_raise(label_i + " = " + self.label + ".frames[" + str(index) + "]")
                remoteframe = RemoteFrame(label_i, self.session)

                # Get the frame
                remoteframes.append(remoteframe)

        # Calculate convolved frames (if necessary)
        if len(for_convolution) > 0:
            # Debugging
            log.debug(str(len(for_convolution)) + " filters require spectral convolution")
            convolved_frames = self.convolve_with_filters(for_convolution, nprocesses=nprocesses, check_previous_sessions=check_previous_sessions)
        else:
            # Debugging
            log.debug("Spectral convolution will be used for none of the filters")
            convolved_frames = []

        # Add the convolved frames
        for fltr, remoteframe in zip(for_convolution, convolved_frames):

            # Set the remote frame
            index = filters.index(fltr)
            remoteframes[index] = remoteframe

        # Return the list of remote frames
        return remoteframes

    # -----------------------------------------------------------------

    def convolve_with_filters(self, filters, nprocesses=8, check_previous_sessions=False):

        """
        This function ...
        :param filters:
        :param nprocesses:
        :param check_previous_sessions:
        :return:
        """

        # Initialize filter list remotely
        self.session.send_line_and_raise("filters = []")

        # Reconstruct the list of filters remotely
        for fltr in filters: self.session.send_line_and_raise("filters.append(BroadBandFilter('" + str(fltr) + "'))")

        # Initialize a list with remoteframes
        remoteframes = []

        # Do the convolution remotely
        self.session.send_line_and_raise("filterconvolvedframes = " + self.label + ".convolve_with_filters(filters, nprocesses=" + str(nprocesses) + ", check_previous_sessions=" + str(check_previous_sessions) + ")", timeout=None, show_output=True)

        # Create a remoteframe pointing to each of the frames in 'filterconvolvedframes'
        last_label = None
        for i in range(len(filters)):

            # Assign a remote label to this result frame
            label_i = get_new_label("Frame", self.session)

            # Check
            if label_i == last_label: raise RuntimeError("Something went wrong: previous frame was not created")

            # Do the assignment remotely
            self.session.send_line_and_raise(label_i + " = filterconvolvedframes[" + str(i) + "]")

            # Create remoteframe and add it to the list
            remoteframe = RemoteFrame(label_i, self.session)
            remoteframes.append(remoteframe)

            # Set last label
            last_label = label_i

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
        self.session.send_line_and_raise(self.label + ".to_wavelength_density('" + new_unit + "', '" + wavelength_unit + "')")

# -----------------------------------------------------------------
