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
import re
from subprocess import call, check_output

# Import the relevant PTS classes and modules
from ...magic.core.frame import Frame
from ..basics.remote import Remote
from . import filesystem, time
from .logging import log

# -----------------------------------------------------------------

# From: http://apple.stackexchange.com/questions/128297/how-to-create-a-vpn-connection-via-terminal
def get_vpn_services():

    """
    This function ...
    :return:
    """

    vpns_string = check_output(["scutil", "--nc", "list"]) # lists all VPN services
    vpns = re.findall('"(.+)"', vpns_string) # service names are double-quoted

    return vpns

# -----------------------------------------------------------------

# From: http://apple.stackexchange.com/questions/128297/how-to-create-a-vpn-connection-via-terminal
def connect_to_vpn(service, user_name=None, password=None, secret=None):

    """
    This function ...
    :param service:
    :param user_name:
    :param password:
    :param secret:
    :return:
    """

    command = ["scutil", "--nc", "start", service]

    if user_name is not None: command += ["--user", user_name]
    if password is not None: command += ["--password", password]
    if secret is not None: command += ["--secret", secret]

    call(command)

# -----------------------------------------------------------------

def remote_convolution(image, kernel_path, kernel_fwhm, host_id):

    """
    This function ...
    :param image:
    :param kernel_path:
    :param kernel_fwhm:
    :param host_id:
    """

    # Debugging
    log.debug("Logging in to remote host ...")

    # Create a remote instance for the specified host ID
    remote = Remote()
    remote.setup(host_id)

    # Debugging
    log.debug("Creating temporary directory remotely ...")

    # Create a temporary directory to do the convolution
    remote_home_directory = remote.home_directory
    remote_temp_path = filesystem.join(remote_home_directory, time.unique_name("convolution"))
    remote.create_directory(remote_temp_path)

    # Debugging
    log.debug("Uploading the kernel to the remote directory ...")

    # Upload the kernel FITS file to the remote directory
    remote_kernel_path = filesystem.join(remote_temp_path, "kernel.fits")
    remote.upload(kernel_path, remote_temp_path, new_name="kernel.fits", compress=True, show_output=True)

    # Debugging
    log.debug("Uploading the image frames to the remote directory ...")

    # Upload the image to the remote directory
    #local_image_path = self.full_output_path("converted_unit.fits")
    #remote_image_path = filesystem.join(remote_temp_path, "image.fits")
    #remote.upload(local_image_path, remote_temp_path, new_name="image.fits", compress=True)

    # Create a temporary directory locally to contain the frames
    local_temp_path = filesystem.join(filesystem.home(), time.unique_name("convolution"))
    filesystem.create_directory(local_temp_path)

    # Save the frames
    local_frame_paths = []
    for frame_name in image.frames:
        frame_path = filesystem.join(local_temp_path, frame_name + ".fits")
        image.frames[frame_name].save(frame_path)
        local_frame_paths.append(frame_path)

    # Upload the frames
    remote_frame_paths = []
    for local_frame_path in local_frame_paths:

        # Determine the name of the local frame file
        frame_file_name = filesystem.name(local_frame_path)

        # Debugging
        log.debug("Uploading the " + filesystem.strip_extension(frame_file_name) + " frame ...")

        # Upload the frame file
        remote_frame_path = filesystem.join(remote_temp_path, frame_file_name)
        remote.upload(local_frame_path, remote_temp_path, new_name=frame_file_name, compress=True, show_output=True)
        remote_frame_paths.append(remote_frame_path)

    # Debugging
    log.debug("Creating a python script to perform the convolution remotely ...")

    # Create a python script that does the convolution
    #script_file = tempfile.NamedTemporaryFile()
    #local_script_path = script_file.name

    local_script_path = filesystem.join(local_temp_path, "convolve.py")
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
    script_file.write("from pts.core.tools.logging import log\n")
    script_file.write("\n")
    script_file.write("# Inform the user\n")
    script_file.write("log.info('Opening the kernel frame ...')\n")
    script_file.write("\n")
    script_file.write("# Open the kernel frame\n")
    script_file.write("kernel = Frame.from_file('" + remote_kernel_path + "')\n")
    script_file.write("\n")
    script_file.write("# Set the FWHM of the kernel\n")
    script_file.write("fwhm = " + str(kernel_fwhm.to("arcsec").value) + " * Unit('arcsec')\n")
    script_file.write("kernel.fwhm = fwhm\n")
    script_file.write("\n")
    #script_file.write("# Open the image\n")
    #script_file.write("image = Image.from_file(" + remote_image_path + ")\n")
    #script_file.write("\n")
    #script_file.write("# Do the convolution")
    #script_file.write("image.convolve(kernel)\n")
    #script_file.write("\n")
    for remote_frame_path in remote_frame_paths:

        frame_name = filesystem.strip_extension(filesystem.name(remote_frame_path))

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
        script_file.write("convolved = frame.convolved(kernel)\n")
        script_file.write("convolved.save('" + remote_frame_path + "')\n") # overwrite the frame
        script_file.write("\n")

        #script_file.write("# Save the image\n")
        #script_file.write("image.save(" + remote_image_path + ")\n")

    # Write to disk
    #script_file.flush()
    script_file.close()

    # Debugging
    log.debug("Uploading the python script ...")

    # Upload the script file
    remote_script_path = filesystem.join(remote_temp_path, "convolve.py")
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
    #remote.download(remote_image_path, filesystem.directory_of(local_result_path), new_name="convolved.fits", compress=True)

    for remote_frame_path in remote_frame_paths:

        # Determine the name of the local frame file
        frame_file_name = filesystem.name(remote_frame_path)

        # Debugging
        log.debug("Downloading the " + filesystem.strip_extension(frame_file_name) + " frame ...")

        # Download
        remote.download(remote_frame_path, local_temp_path, new_name=frame_file_name, compress=True, show_output=True)

    # Remove the temporary directory on the remote's filesystem
    remote.remove_directory(remote_temp_path)

    # Load the result
    #self.image = Image.from_file(local_result_path)

    for frame_name in image.frames.keys():
        local_frame_path = filesystem.join(local_temp_path, frame_name + ".fits")
        image.frames[frame_name] = Frame.from_file(local_frame_path)

    # Remove the local temporary directory
    filesystem.remove_directory(local_temp_path)

# -----------------------------------------------------------------

def remote_convolution2(frame, kernel, host_id):

    """
    This function ...
    :param frame:
    :param kernel:
    :param host_id:
    :return:
    """

    # Debugging
    log.debug("Logging in to remote host ...")

    # Create a remote instance for the specified host ID
    remote = Remote()
    remote.setup(host_id)

    # Debugging
    log.debug("Creating temporary directory remotely ...")

    # Create a temporary directory to do the convolution
    remote_home_directory = remote.home_directory
    remote_temp_path = filesystem.join(remote_home_directory, time.unique_name("convolution"))
    remote.create_directory(remote_temp_path)

    # Debugging
    log.debug("Creating local temporary directory ...")

    # Create a temporary directory locally to contain the frames
    local_temp_path = filesystem.join(filesystem.home(), time.unique_name("convolution"))
    filesystem.create_directory(local_temp_path)

    # Debugging
    log.debug("Writing the frame to the temporary directory ...")

    # Write the frame
    local_frame_path = filesystem.join(local_temp_path, frame.name + ".fits")
    frame.save(local_frame_path)

    # Debugging
    log.debug("Writing the kernel to the temporary directory ...")

    # Write the kernel
    local_kernel_path = filesystem.join(local_temp_path, "kernel.fits")
    kernel.save(local_kernel_path)

    # Debugging
    log.debug("Uploading the frame to the remote directory ...")

    # Upload the frame file
    remote_frame_path = filesystem.join(remote_temp_path, frame.name)
    remote.upload(local_frame_path, remote_temp_path, new_name=frame.name, compress=True, show_output=True)

    # Debugging
    log.debug("Uploading the kernel to the remote directory ...")

    # Upload the kernel FITS file to the remote directory
    remote_kernel_path = filesystem.join(remote_temp_path, "kernel.fits")
    remote.upload(local_kernel_path, remote_temp_path, new_name="kernel.fits", compress=True, show_output=True)

    # Debugging
    log.debug("Creating a python script to perform the convolution remotely ...")

    # Create the script
    local_script_path = filesystem.join(local_temp_path, "convolve.py")
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
    script_file.write("convolved = frame.convolved(kernel)\n")
    script_file.write("convolved.save('" + remote_frame_path + "')\n") # overwrite the frame

    # Write to disk
    script_file.close()

    # Debugging
    log.debug("Uploading the python script ...")

    # Upload the script file
    remote_script_path = filesystem.join(remote_temp_path, "convolve.py")
    remote.upload(local_script_path, remote_temp_path, new_name="convolve.py", show_output=True)

    # Debugging
    log.debug("Executing the script remotely ...")

    # Execute the script file remotely
    remote.execute("python " + remote_script_path, output=False, show_output=True)

    # Debugging
    log.debug("Downloading the result ...")

    # Determine the name of the local frame file
    frame_file_name = filesystem.name(remote_frame_path)

    # Debugging
    log.debug("Downloading the " + filesystem.strip_extension(frame_file_name) + " frame ...")

    # Download
    remote.download(remote_frame_path, local_temp_path, new_name=frame_file_name, compress=True, show_output=True)

    # Remove the temporary directory on the remote's filesystem
    remote.remove_directory(remote_temp_path)

    # Load the convolved frame
    convolved = Frame.from_file(local_frame_path)

    # Remove the local temporary directory
    filesystem.remove_directory(local_temp_path)

    # Return the convolved frame
    return convolved

# -----------------------------------------------------------------
