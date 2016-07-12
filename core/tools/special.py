#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.special Special functions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...magic.core.frame import Frame
from ..basics.remote import Remote, connected_remotes
from . import time
from . import filesystem as fs
from .logging import log

# -----------------------------------------------------------------

def remote_convolution(image, kernel, host_id):

    """
    This function ...
    :param image:
    :param kernel:
    :param host_id:
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

    # Debugging
    log.debug("Creating temporary directory remotely ...")

    # Create a temporary directory to do the convolution
    remote_home_directory = remote.home_directory
    remote_temp_path = fs.join(remote_home_directory, time.unique_name("convolution"))
    remote.create_directory(remote_temp_path)

    # Debugging
    #log.debug("Uploading the kernel to the remote directory ...")

    # Upload the kernel FITS file to the remote directory
    #remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
    #remote.upload(kernel_path, remote_temp_path, new_name="kernel.fits", compress=True, show_output=True)

    # Debugging
    log.debug("Creating a local temporary directory ...")

    # Create a temporary directory locally to contain the frames
    local_temp_path = fs.join(fs.home(), time.unique_name("convolution"))
    fs.create_directory(local_temp_path)

    # Debugging
    log.debug("Saving the image frames to the temporary directory ...")

    # Save the frames
    local_frame_paths = []
    constant_frames = []
    for frame_name in image.frames:
        frame_path = fs.join(local_temp_path, frame_name + ".fits")

        # Only upload and convolve non-constant frames
        if not image.frames[frame_name].is_constant():

            image.frames[frame_name].save(frame_path)
            local_frame_paths.append(frame_path)

        else:

            log.debug("The " + frame_name + " frame is constant, so this won't be uploaded and convolved")
            constant_frames.append(frame_name)

    # Debugging
    log.debug("Saving the kernel to the temporary directory ...")
    local_kernel_path = fs.join(local_temp_path, "kernel.fits")
    kernel.save(local_kernel_path)

    # Debugging
    log.debug("Uploading the image frames to the remote directory ...")

    # Upload the frames
    remote_frame_paths = []
    for local_frame_path in local_frame_paths:

        # Determine the name of the local frame file
        frame_file_name = fs.name(local_frame_path)

        # Debugging
        log.debug("Uploading the " + fs.strip_extension(frame_file_name) + " frame ...")

        # Upload the frame file
        remote_frame_path = fs.join(remote_temp_path, frame_file_name)
        remote.upload(local_frame_path, remote_temp_path, new_name=frame_file_name, compress=True, show_output=True)
        remote_frame_paths.append(remote_frame_path)

    # Debugging
    log.debug("Uploading the kernel to the remote directory ...")

    # Upload the kernel
    remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
    remote.upload(local_kernel_path, remote_temp_path, new_name="kernel.fits", compress=True, show_output=True)

    # Debugging
    log.debug("Creating a python script to perform the convolution remotely ...")

    # Create a python script that does the convolution
    #script_file = tempfile.NamedTemporaryFile()
    #local_script_path = script_file.name

    local_script_path = fs.join(local_temp_path, "convolve.py")
    script_file = open(local_script_path, 'w')

    script_file.write("#!/usr/bin/env python\n")
    script_file.write("# -*- coding: utf8 -*-\n")
    script_file.write("\n")
    script_file.write("# Import astronomical modules\n")
    script_file.write("from astropy.units import Unit\n")
    script_file.write("\n")
    script_file.write("# Import the relevant PTS classes and modules\n")
    script_file.write("from pts.magic.core.frame import Frame\n")
    script_file.write("from pts.magic.core.image import Image\n")
    script_file.write("from pts.magic.core.kernel import ConvolutionKernel\n")
    script_file.write("from pts.core.tools.logging import log\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Opening the kernel frame ...')\n")
    script_file.write("\n")
    script_file.write("# Open the kernel\n")
    script_file.write("kernel = ConvolutionKernel.from_file('" + remote_kernel_path + "')\n")
    script_file.write("\n")
    for remote_frame_path in remote_frame_paths:

        frame_name = fs.strip_extension(fs.name(remote_frame_path))

        script_file.write("# Inform the user\n")
        script_file.write("log.info('Opening the " + frame_name + " frame ...')\n")
        script_file.write("\n")
        script_file.write("# Open the frame\n")
        script_file.write("frame = Frame.from_file('" + remote_frame_path + "')\n")
        script_file.write("\n")
        script_file.write("# Inform the user\n")
        script_file.write("log.info('Convolving the " + frame_name + " frame ...')\n")
        script_file.write("\n")
        script_file.write("# Do the convolution and save the result\n")
        script_file.write("frame.convolve(kernel, allow_huge=True)\n")
        script_file.write("frame.save('" + remote_frame_path + "')\n") # overwrite the frame
        script_file.write("\n")

        #script_file.write("# Save the image\n")
        #script_file.write("image.save(" + remote_image_path + ")\n")

    # Write to disk
    #script_file.flush()
    script_file.close()

    # Debugging
    log.debug("Uploading the python script ...")

    # Upload the script file
    remote_script_path = fs.join(remote_temp_path, "convolve.py")
    remote.upload(local_script_path, remote_temp_path, new_name="convolve.py", show_output=True)

    # Close the local script (it is automatically removed)
    #script_file.close()

    # Debugging
    log.debug("Executing the script remotely ...")

    # Execute the script file remotely
    remote.execute("python " + remote_script_path, output=False, show_output=True)

    # Debugging
    log.debug("Downloading the results ...")

    # Download the resulting FITS file (the convolved image)
    #local_result_path = self.full_output_path("convolved.fits")
    #remote.download(remote_image_path, fs.directory_of(local_result_path), new_name="convolved.fits", compress=True)

    for remote_frame_path in remote_frame_paths:

        # Determine the name of the local frame file
        frame_file_name = fs.name(remote_frame_path)

        # Debugging
        log.debug("Downloading the " + fs.strip_extension(frame_file_name) + " frame ...")

        # Download
        remote.download(remote_frame_path, local_temp_path, new_name=frame_file_name, compress=True, show_output=True)

    # Remove the temporary directory on the remote's filesystem
    remote.remove_directory(remote_temp_path)

    # Load the result
    #self.image = Image.from_file(local_result_path)

    for frame_name in image.frames.keys():
        if frame_name in constant_frames: continue # Skip constant frames, these are not convolved
        local_frame_path = fs.join(local_temp_path, frame_name + ".fits")
        image.frames[frame_name] = Frame.from_file(local_frame_path)

    # Remove the local temporary directory
    fs.remove_directory(local_temp_path)

# -----------------------------------------------------------------

def remote_convolution_frame(frame, kernel_path, host_id):

    """
    This function ...
    :param frame:
    :param kernel_path:
    :param host_id:
    :return:
    """

    # Check whether the frame is constant. If it is, we don't have to convolve!
    if frame.is_constant(): return frame.copy()

    # Check whether we are already connected to the specified remote host
    if host_id in connected_remotes and connected_remotes[host_id] is not None:
        remote = connected_remotes[host_id]
    else:

        # Debugging
        log.debug("Logging in to remote host ...")

        # Create a remote instance for the specified host ID
        remote = Remote()
        remote.setup(host_id)

    # Debugging
    log.debug("Creating temporary directory remotely ...")

    # Create a temporary directory to do the convolution
    remote_home_directory = remote.home_directory
    remote_temp_path = fs.join(remote_home_directory, time.unique_name("convolution"))
    remote.create_directory(remote_temp_path)

    # Debugging
    log.debug("Creating local temporary directory ...")

    # Create a temporary directory locally to contain the frames
    local_temp_path = fs.join(fs.home(), time.unique_name("convolution"))
    fs.create_directory(local_temp_path)

    # Debugging
    log.debug("Writing the frame to the temporary directory ...")

    # Write the frame
    local_frame_path = fs.join(local_temp_path, frame.name + ".fits")
    frame.save(local_frame_path)

    # Debugging
    #log.debug("Writing the kernel to the temporary directory ...")

    # Write the kernel
    #local_kernel_path = fs.join(local_temp_path, "kernel.fits")
    #kernel.save(local_kernel_path)

    # Debugging
    log.debug("Uploading the frame to the remote directory ...")

    # Upload the frame file
    remote_frame_path = fs.join(remote_temp_path, frame.name)
    remote.upload(local_frame_path, remote_temp_path, new_name=frame.name, compress=True, show_output=True)

    # Debugging
    #log.debug("Uploading the kernel to the remote directory ...")

    # Upload the kernel FITS file to the remote directory
    #remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
    #remote.upload(local_kernel_path, remote_temp_path, new_name="kernel.fits", compress=True, show_output=True)

    # Debugging
    log.debug("Uploading the kernel to the remote directory ...")

    # Upload the kernel FITS file to the remote directory
    remote_kernel_path = fs.join(remote_temp_path, "kernel.fits")
    remote.upload(kernel_path, remote_temp_path, new_name="kernel.fits", compress=True, show_output=True)

    # Debugging
    log.debug("Creating a python script to perform the convolution remotely ...")

    # Create the script
    local_script_path = fs.join(local_temp_path, "convolve.py")
    script_file = open(local_script_path, 'w')

    script_file.write("#!/usr/bin/env python\n")
    script_file.write("# -*- coding: utf8 -*-\n")
    script_file.write("\n")
    script_file.write("# Import the relevant PTS classes and modules\n")
    script_file.write("from pts.magic.core.frame import Frame\n")
    script_file.write("from pts.core.tools.logging import log\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Opening the kernel frame ...')\n")
    script_file.write("\n")
    script_file.write("# Open the kernel frame\n")
    script_file.write("kernel = Frame.from_file('" + remote_kernel_path + "')\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Opening the frame ...')\n")
    script_file.write("\n")
    script_file.write("# Open the frame\n")
    script_file.write("frame = Frame.from_file('" + remote_frame_path + "')\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Convolving the frame ...')\n")
    script_file.write("\n")
    script_file.write("# Do the convolution and save the result\n")
    script_file.write("convolved = frame.convolved(kernel, allow_huge=True)\n")
    script_file.write("convolved.save('" + remote_frame_path + "')\n") # overwrite the frame

    # Write to disk
    script_file.close()

    # Debugging
    log.debug("Uploading the python script ...")

    # Upload the script file
    remote_script_path = fs.join(remote_temp_path, "convolve.py")
    remote.upload(local_script_path, remote_temp_path, new_name="convolve.py", show_output=True)

    # Debugging
    log.debug("Executing the script remotely ...")

    # Execute the script file remotely
    remote.execute("python " + remote_script_path, output=False, show_output=True)

    # Debugging
    log.debug("Downloading the result ...")

    # Determine the name of the local frame file
    frame_file_name = fs.name(remote_frame_path)

    # Debugging
    log.debug("Downloading the " + fs.strip_extension(frame_file_name) + " frame ...")

    # Download
    remote.download(remote_frame_path, local_temp_path, new_name=frame_file_name, compress=True, show_output=True)

    # Remove the temporary directory on the remote's filesystem
    remote.remove_directory(remote_temp_path)

    # Load the convolved frame
    convolved = Frame.from_file(local_frame_path)

    # Remove the local temporary directory
    fs.remove_directory(local_temp_path)

    # Return the convolved frame
    return convolved

# -----------------------------------------------------------------

def remote_filter_convolution_no_pts(host_id, datacube_path, wavelengths, filters):

    """
    This function ...
    :param host_id:
    :param datacube_path:
    :param wavelengths:
    :param filters:
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


    # Debugging
    log.debug("Creating temporary directory remotely ...")

    # Create a temporary directory to do the convolution
    remote_home_directory = remote.home_directory
    remote_temp_path = fs.join(remote_home_directory, time.unique_name("filter-convolution"))
    remote.create_directory(remote_temp_path)


    # Debugging
    log.debug("Creating local temporary directory ...")

    # Create a temporary directory locally to contain the frames
    local_temp_path = fs.join(fs.home(), time.unique_name("filter-convolution"))
    fs.create_directory(local_temp_path)

    integrated_transmissions = dict()

    # Loop over the filters
    for fltr in filters:

        # Get the transmission data
        fltr_wavelengths = fltr._Wavelengths
        fltr_transmission = fltr._Transmission
        fltr_integrated_transmission = fltr._IntegratedTransmission
        integrated_transmissions[fltr.name] = fltr_integrated_transmission

        # Save the transmission data
        path = fs.join(local_temp_path, "transmission__" + str(fltr) + ".dat")
        np.savetxt(path, (fltr_wavelengths, fltr_transmission))

    #print(integrated_transmissions)
    #print(local_temp_path)

    integrated_path = fs.join(local_temp_path, "integrated_transmissions.txt")
    with open(integrated_path, 'w') as integrated_trans_file:
        for fltr_name in integrated_transmissions:
            integrated_trans_file.write(fltr_name + ": " + str(integrated_transmissions[fltr_name]) + "\n")

    # NOT FINISHED ...

# -----------------------------------------------------------------

def remote_filter_convolution(host_id, datacube_path, wavelengths, filters, keep_output=False):

    """
    This function ...
    :param host_id:
    :param datacube_path:
    :param wavelengths:
    :param filters:
    :param keep_output:
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

    # Debugging
    log.debug("Creating temporary directory remotely ...")

    # Create a temporary directory to do the convolution
    remote_home_directory = remote.home_directory
    remote_temp_path = fs.join(remote_home_directory, time.unique_name("filter-convolution"))
    remote.create_directory(remote_temp_path)

    # Debugging
    log.debug("Creating local temporary directory ...")

    # Create a temporary directory locally to contain the frames
    local_temp_path = fs.join(fs.home(), time.unique_name("filter-convolution"))
    fs.create_directory(local_temp_path)

    # Debugging
    log.debug("Uploading the datacube to the temporary remote directory ...")

    # Upload the frame file
    datacube_name = fs.name(datacube_path)
    remote_datacube_path = fs.join(remote_temp_path, datacube_name)
    remote.upload(datacube_path, remote_temp_path, compress=True, show_output=True)

    # Debugging
    log.debug("Writing the wavelengths to the temporary local directory ...")
    local_wavelengths_path = fs.join(local_temp_path, "wavelengths.txt")
    np.savetxt(local_wavelengths_path, wavelengths)

    # Debugging
    log.debug("Uploading the wavelengths file to the remote directory ...")

    # Upload the kernel FITS file to the remote directory
    remote_wavelengths_path = fs.join(remote_temp_path, "wavelengths.txt")
    remote.upload(local_wavelengths_path, remote_temp_path, compress=True, show_output=True)

    # Debugging
    log.debug("Creating a python script to perform the filter convolution remotely ...")

    # Create the script
    local_script_path = fs.join(local_temp_path, "make_images.py")
    script_file = open(local_script_path, 'w')

    script_file.write("#!/usr/bin/env python\n")
    script_file.write("# -*- coding: utf8 -*-\n")
    script_file.write("\n")
    script_file.write("# Import standard modules\n")
    script_file.write("import numpy as np\n")
    script_file.write("\n")
    script_file.write("# Import the relevant PTS classes and modules\n")
    script_file.write("from pts.magic.core.image import Image\n")
    script_file.write("from pts.magic.core.frame import Frame\n")
    script_file.write("from pts.core.basics.filter import Filter\n")
    script_file.write("from pts.core.tools.logging import log\n")
    script_file.write("from pts.core.tools import filesystem as fs\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Loading the datacube ...')\n")
    script_file.write("\n")
    script_file.write("# Open the datacube as an Image\n")
    script_file.write("datacube = Image.from_file('" + remote_datacube_path + "', always_call_first_primary=False)\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Loading the wavelengths ...')\n")
    script_file.write("\n")
    script_file.write("# Load the wavelengths from the text file\n")
    script_file.write("wavelengths = np.loadtxt('" + remote_wavelengths_path + "')\n")
    script_file.write("\n")
    script_file.write("# Convert the frames from neutral surface brightness to wavelength surface brightness\n")
    script_file.write("for l in range(len(wavelengths)):\n")
    script_file.write("\n")
    script_file.write("    # Get the wavelength\n")
    script_file.write("    wavelength = wavelengths[l]\n")
    script_file.write("\n")
    script_file.write("    # Determine the name of the frame in the datacube\n")
    script_file.write("    frame_name = 'frame' + str(l)\n")
    script_file.write("\n")
    script_file.write("    # Divide this frame by the wavelength in micron\n")
    script_file.write("    datacube.frames[frame_name] /= wavelength\n")
    script_file.write("\n")
    script_file.write("    # Set the new unit\n")
    script_file.write("    datacube.frames[frame_name].unit = 'W / (m2 * arcsec2 * micron)'\n")
    script_file.write("\n")
    script_file.write("# Convert the datacube to a numpy array where wavelength is the third dimension\n")
    script_file.write("fluxdensities = datacube.asarray()\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Creating the filters ...')\n")
    script_file.write("\n")
    script_file.write("filters = dict()\n")
    script_file.write("\n")
    for filter_name in filters:
        fltr = filters[filter_name]
        script_file.write("# Inform the user\n")
        script_file.write("log.info('Creating the " + str(fltr) + " filter')\n")
        script_file.write("\n")
        script_file.write("fltr = Filter.from_string('" + str(fltr) + "')\n")
        script_file.write("filters['" + filter_name + "'] = fltr\n")
        script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Performing the filter convolutions ...')\n")
    script_file.write("\n")
    script_file.write("# Loop over the filters, perform the convolution\n")
    script_file.write("for filter_name in filters:\n")
    script_file.write("\n")
    script_file.write("    log.info('Making the observed image for the ' + str(fltr) + ' filter ...')\n")
    script_file.write("    fltr = filters[filter_name]\n")
    script_file.write("    data = fltr.convolve(wavelengths, fluxdensities)\n")
    script_file.write("    frame = Frame(data)\n")
    script_file.write("    frame.unit = 'W/(m2 * arcsec2 * micron)'\n")
    script_file.write("    path = fs.join('" + remote_temp_path + "', filter_name + '.fits')\n")
    script_file.write("    frame.save(path)\n")

    # Write to disk
    script_file.close()

    # Debugging
    log.debug("Uploading the python script ...")

    # Upload the script file
    remote_script_path = fs.join(remote_temp_path, "make_images.py")
    remote.upload(local_script_path, remote_temp_path, new_name="make_images.py", show_output=True)

    # Debugging
    log.debug("Executing the script remotely ...")

    # Execute the script file remotely
    remote.execute("python " + remote_script_path, output=False, show_output=True)

    # Remove the datacube in the remote directory
    remote.remove_file(remote_datacube_path)

    # Debugging
    log.debug("Downloading the convolved frames ...")

    # Download
    local_downloaded_temp_path = fs.join(fs.home(), fs.name(remote_temp_path))
    fs.create_directory(local_downloaded_temp_path)
    remote.download(remote_temp_path, local_downloaded_temp_path, compress=True, show_output=True)

    # Remove the temporary directory on the remote's filesystem
    remote.remove_directory(remote_temp_path)

    # Remove the local temporary directory
    fs.remove_directory(local_temp_path)

    # Create a dictionary to contain the frames
    frames = dict()

    # Loop over the filters, load the frame
    for filter_name in filters:

        # Determine the path to the resulting FITS file
        path = fs.join(local_downloaded_temp_path, filter_name + ".fits")

        # Check whether the frame exists
        if not fs.is_file(path): raise RuntimeError("The image for filter " + str(filters[filter_name]) + " is missing")

        # Load the FITS file
        frame = Frame.from_file(path)

        # Add the frame to the dictionary
        frames[filter_name] = frame

    # Remove the downloaded temporary directory
    if not keep_output: fs.remove_directory(local_downloaded_temp_path)

    # Return the dictionary of frames
    return frames

# -----------------------------------------------------------------
