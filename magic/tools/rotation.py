from __future__ import print_function

# Import standard modules
import os
import numpy as np
from scipy import ndimage
from scipy import misc

# Import astronomical modules
from astropy.io import fits

# This script rotates a FITS image around the center of the frame over an arbitrary angle.
# It updates the header accordingly.
# Written by Sebastien Viaene, Universiteit Gent, May 13 2015

def main():

    # Get all the images in the ScanMode directory
    in_directory = os.path.join("/home/mbaes/DustPedia", "ScanMode")
    out_directory = os.path.join(os.getcwd(), "PACS", "ScanMode-rotated")
    
    # Loop over all images
    for name in os.listdir(in_directory):

        # Skip non-FITS files
        if not name.endswith(".fits") or name.startswith("."): continue

        # Determine the full path to the file
        in_path = os.path.join(in_directory, name)
        out_path = os.path.join(out_directory, name)
    
        # Inform the user
        print("Processing ", in_path)
    
        # Load in image
        hdu = fits.open(in_path)[0]

        # Get the rotation angle
        angle = hdu.header["CROTA2"]

        # Rotate the header
        new_header = rotate_header(hdu.header, angle)

        print("Shape:", hdu.shape)

        # Get the number of frames
        n_frames = hdu.shape[0] if len(hdu.shape) == 3 else 1
        
        # Initialize a new image
        new_image = [None] * n_frames

        # Only one frame
        if n_frames == 1: new_image = rotate_frame(hdu.data[0:,0:], angle)
        
        # If there are multiple frames
        else:

            # Rotate each of the image frames
            for i in range(n_frames):
            
                # Create a rotated version of the frame and add it to the new image
                new_frame = rotate_frame(hdu.data[i][0:,0:], angle)
                new_image[i] = new_frame

        # Write out rotated image.
        rot_hdu = fits.PrimaryHDU(new_image, new_header)
        rot_hdu.writeto(out_path, clobber=True)

def rotate_header(header, angle):
    
    """
    This function rotates the header
    """
    
    new_header = header
    
    # Check if a rotation matrix element exists
    matrix = True
    try:
        cd1_1 = np.float(header["CD1_1"])
    except:
        matrix = False

    theta = angle * np.pi / 180.
    rot_matrix = np.array( [ [ np.cos(theta), np.sin(theta)],
                            [-1.*np.sin(theta), np.cos(theta)] ] )

    center = np.array([(header['NAXIS1'] - 1)/2., (header['NAXIS2'] - 1)/2. ])
    
    try:
        crpix = np.array([header['CRPIX1'], header['CRPIX2']])
    except:
        crpix = center
        header.append( fits.Card('CRPIX1', crpix[0], 'Reference pixel on this axis'), end=True)
        header.append( fits.Card('CRPIX2', crpix[1], 'Reference pixel on this axis'), end=True)

    ncrpix = (crpix-1-center).dot(rot_matrix.T) + 1
    ncrpix += center

    new_header["CRPIX1"] = ncrpix[0]
    new_header["CRPIX2"] = ncrpix[1]

    if matrix:

        try:
            cd1_2 = np.float(header["CD1_2"])
        except:
            cd1_2 = 0.
            header.append(fits.Card('CD1_2', cd1_2, 'Rotation matrix element 1_2'), end=True)

        try:
            cd2_1 = np.float(header["CD2_1"])
        except:
            cd2_1 = 0.
            header.append(fits.Card('CD2_1', cd2_1, 'Rotation matrix element 2_1'), end=True)

        try:
            cd2_2 = np.float(header["CD2_2"])
        except:
            cd2_2 = 0.
            header.append(fits.Card('CD2_2', cd2_2, 'Rotation matrix element 2_2'), end=True)

        cd = np.array([[cd1_1, cd1_2], [cd2_1, cd2_2]])

        newcd = rot_matrix.dot(cd)
        new_header["CD1_1"] = newcd[0,0]
        new_header["CD1_2"] = newcd[0,1]
        new_header["CD2_1"] = newcd[1,0]
        new_header["CD2_2"] = newcd[1,1]
    
    else:
        
        #try:
        #    new_header["CROTA1"] = -1.*angle
        #except:
        #    new_header.append(fits.Card('CROTA1', -1.*angle, 'Rotation parameter'), end=True)

        #try:
        #    new_header["CROTA2"] = -1.*angle
        #except:
        #   new_header.append( fits.Card('CROTA2', -1.*angle, 'Rotation parameter'), end=True)

        new_header["CROTA2"] = 0.0

    return new_header

def rotate_frame(frame, angle):
    
    # Perform the image rotation and update the fits header
    #frame[np.isnan(frame)] = 0.0
    new_frame = ndimage.interpolation.rotate(frame, angle, reshape=False, order=1, mode='constant', cval=float('nan'))
    
    #new_frame = misc.imrotate(frame, angle, interp="bilinear")
    
    # Return the rotated frame
    return new_frame

def setRelativeWCS(hdr, xcenter, ycenter, pxlScale):

    # remove anything but the basic header

    hdr.append(fits.Card('BOTTOM', 'END', 'End of header'), end=True)
    i= 5
    while hdr[i] != 'END':
        del hdr[i]
    del hdr['BOTTOM']

    # add relative coordinate info
    hdr.append(fits.Card('CRPIX1', xcenter, 'Reference pixel on this axis'), end=True)
    hdr.append(fits.Card('CDELT1', pxlScale, 'Coordinate increment along  this axis'), end=True)
    hdr.append(fits.Card('CRVAL1', 0., 'World coordinate on this axis'), end=True)
    hdr.append(fits.Card('CTYPE1', 'LINEAR', 'WCS projection type for this axis'), end=True)
    hdr.append(fits.Card('CRPIX2', ycenter, 'Reference pixel on this axis'), end=True)
    hdr.append(fits.Card('CDELT2', pxlScale, 'Coordinate increment along  this axis'), end=True)
    hdr.append(fits.Card('CRVAL2', 0., 'World coordinate on this axis'), end=True)
    hdr.append(fits.Card('CTYPE2', 'LINEAR', 'WCS projection type for this axis'), end=True)

    return hdr

def cropImage(im, hdr, region):
    
    xmin, ymin, xmax, ymax = readBox(region)
    cim = im[ymin-1:ymax-1,xmin-1:xmax-1]
    hdr['CRPIX1'] = hdr['CRPIX1'] - xmin + 1
    hdr['CRPIX2'] = hdr['CRPIX2'] - ymin + 1
    return cim, hdr

def readBox(region):
    
    with open(region) as f:
        for _ in xrange(3):
            line = f.readline()
        line = f.readline().replace('(',',').replace(')',',').split(',')
    
    xc = float(line[1])
    yc = float(line[2])
    xsize = float(line[3])
    ysize = float(line[4])

    return int(round(xc-xsize/2.)), int(round(yc-ysize/2.)), int(round(xc+xsize/2.)), int(round(yc+ysize/2.))
