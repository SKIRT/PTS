#!/bin/bash
# (use "chmod +rx scriptname" to make script executable)
#
# This script publishes the PTS documentation to the public web host.
#
# Instructions:
#   - use on Mac OS X only
#   - before running this script, build the PDF and HTML documentation
#   - you need write access to the public web host
#   - run this script with "git" as default directory

# mount the webhost at a local mount point
mkdir -p ../webhost
mount_smbfs smb://pcamps@webhost.ugent.be/_skirt ../webhost

# copy the downloadable documentation
mkdir -p ../webhost/WWW/downloads
cp -v ../doc/PTS.qch ../webhost/WWW/downloads/

# copy the online html documentation
rm -rf ../webhost/WWW/pts
cp -Rv ../html ../webhost/WWW/pts

# unmount the webhost and remove the mount point
umount ../webhost
rm -d ../webhost
