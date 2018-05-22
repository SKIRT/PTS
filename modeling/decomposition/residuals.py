#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.residuals Contains the DecompositionResidualsCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.modeling.models import Gaussian2D

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.magic.core.source import Source
from pts.magic.tools import statistics, plotting, fitting
from pts.magic.region.ellipse import PixelEllipseRegion
from pts.magic.basics.vector import Extent
from .component import DecompositionComponent
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

class DecompositionResidualsCalculator(DecompositionComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        super(DecompositionResidualsCalculator, self).__init__(*args, **kwargs)

        self.i1_jy = None
        self.disk_jy = None
        self.bulge2d_jy = None
        self.bulge_jy = None
        self.model_jy = None

        # Residuals
        self.bulge2d_residual = None
        self.bulge_residual = None
        self.disk_residual = None
        self.model_residual = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Load the images
        self.load_images()

        # Calculate residuals
        self.calculate_residuals()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        self.load_i1()

        self.load_bulge2d()

        self.load_bulge()

        self.load_disk()

        self.load_model()

    # -----------------------------------------------------------------

    def load_i1(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the IRAC 3.6 micron image ...")

        # Determine the path to the truncated 3.6 micron image
        #path = fs.join(truncation_path, "IRAC I1.fits")

        path = fs.join(self.prep_path, "IRAC I1", "result.fits")
        frame = Frame.from_file(path)

        # Convert the frame to Jy/pix
        conversion_factor = 1.0
        conversion_factor *= 1e6

        # Convert the 3.6 micron image from Jy / sr to Jy / pixel
        pixelscale = frame.average_pixelscale
        pixel_factor = (1.0/pixelscale**2).to("pix2/sr").value
        conversion_factor /= pixel_factor
        frame *= conversion_factor
        frame.unit = "Jy"

        # Set the frame
        self.i1_jy = frame

    # -----------------------------------------------------------------

    def load_bulge2d(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the bulge 2D image ...")

        # Determine the path to the bulge 2D image
        path = fs.join(self.components_images_path, "bulge2D.fits")
        self.bulge2d_jy = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the bulge image ...")

        # Determine the path to the truncated bulge image
        path = fs.join(self.components_images_path, "bulge.fits")
        self.bulge_jy = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the disk image ...")

        # Determine the path to the disk image
        path = fs.join(self.components_images_path, "disk.fits")
        self.disk_jy = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the model image ...")

        # Determine the path to the truncated model image
        path = fs.join(self.components_images_path, "model.fits")
        self.model_jy = Frame.from_file(path)

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Calculate the bulge2D residual
        self.bulge2d_residual = self.i1_jy - self.bulge2d_jy

        # Calculate the bulge residual frame
        self.bulge_residual = self.i1_jy - self.bulge_jy

        # Calculate the disk residual
        self.disk_residual = self.i1_jy - self.disk_jy

        # Calculate the model residual
        self.model_residual = self.i1_jy - self.model_jy

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write ...
        self.write_jy()

        # Write the residual frames
        self.write_residuals()

    # -----------------------------------------------------------------

    def write_jy(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_residuals(self):

        """
        This function ...
        :return:
        """

        # Bulge residual
        path = fs.join(self.components_residuals_path, "bulge.fits")
        self.bulge_residual.saveto(path)

        # Bulge2D residual
        path = fs.join(self.components_residuals_path, "bulge2D.fits")
        self.bulge2d_residual.saveto(path)

        # Disk residual
        path = fs.join(self.components_residuals_path, "disk.fits")
        self.disk_residual.saveto(path)

        # Calculate the model residual frame
        path = fs.join(self.components_residuals_path, "model.fits")
        self.model_residual.saveto(path)

    # -----------------------------------------------------------------

    def old(self):

        """
        This function ...
        :return:
        """

        exit()

        # FWHM of all the images
        fwhm = 11.18 * u("arcsec")
        fwhm_pix = (fwhm / frame.average_pixelscale).to("pix").value
        sigma = fwhm_pix * statistics.fwhm_to_sigma

        # Get the center pixel of the galaxy
        parameters_path = fs.join(components_path, "parameters.dat")
        parameters = load_parameters(parameters_path)
        center = parameters.center.to_pixel(frame.wcs)

        # Create a source around the galaxy center
        ellipse = PixelEllipseRegion(center, 20.0*sigma)
        source = Source.from_ellipse(model_residual, ellipse, 1.5)

        source.estimate_background("polynomial")

        source.plot()

        position = source.center
        model = source.subtracted.fit_model(position, "Gaussian")

        rel_center = center - Extent(source.x_min, source.y_min)
        rel_model = fitting.shifted_model(model, -source.cutout.x_min, -source.cutout.y_min)
        plotting.plot_peak_model(source.cutout, rel_center.x, rel_center.y, rel_model)

        model_fwhm_pix = fitting.fwhm(model)
        model_fwhm = (model_fwhm_pix * frame.average_pixelscale).to("arcsec")

        print("Model FWHM: ", model_fwhm)

        evaluated_model = source.cutout.evaluate_model(model)

        all_residual = Frame(np.copy(model_residual))
        all_residual[source.y_slice, source.x_slice] -= evaluated_model
        all_residual.saveto(fs.join(residuals_path, "all_residual.fits"))

        model = Gaussian2D(amplitude=0.0087509425805, x_mean=center.x, y_mean=center.y, x_stddev=sigma, y_stddev=sigma)
        rel_model = fitting.shifted_model(model, -source.cutout.x_min, -source.cutout.y_min)
        plotting.plot_peak_model(source.cutout, rel_center.x, rel_center.y, rel_model)

        evaluated_model = source.cutout.evaluate_model(model)

        all_residual2 = Frame(np.copy(model_residual))
        all_residual2[source.y_slice, source.x_slice] -= evaluated_model
        all_residual2.saveto(fs.join(residuals_path, "all_residual2.fits"))

# -----------------------------------------------------------------
