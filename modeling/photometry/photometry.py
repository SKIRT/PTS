#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.photometry Contains the PhotoMeter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import PhotometryComponent
from .sedfetching import SEDFetcher
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ..core.sed import ObservedSED
from ...core.basics.errorbar import ErrorBar
from ...core.tools import tables
from ...core.plot.sed import SEDPlotter
from ...magic.misc.kernels import AnianoKernels
from ...magic.core.kernel import ConvolutionKernel

# -----------------------------------------------------------------

# CHRIS:

# For SPIRE: Section 5.2.12 of the SPIRE Handbook (http://herschel.esac.esa.int/Docs/SPIRE/spire_handbook.pdf)
# provides the "official" way of performing aperture corrections; we specifically care about the "Extended source photometry
# (starting from extended source maps)" part. If you're using arbitrarily-large, non-circular apertures, you basically
# have to use maps of the beam profile to work out what fraction of the flux is outside your aperture.
# Those can be retrieved from the SPIRE calibration wiki
# (http://herschel.esac.esa.int/twiki/bin/view/Public/SpirePhotometerBeamProfile2).

# In the meantime, I'm working on sensible automated aperture-corrections - current plan is to assume that the underlying
# flux distribution follows a 2D-sersic profile, and so fit to each source a 2D-sersic convolved with the beam, and
# hence estimate the amount of flux that gets spread beyond the aperture. Hopefully it will be easily applicable
# to your work too.

# For PACS: The most up-to-date document on the PACS calibration wiki
# (http://herschel.esac.esa.int/twiki/pub/Public/PacsCalibrationWeb/bolopsf_22.pdf)
# and its accompanying tar.gz (ftp://ftp.sciops.esa.int/pub/hsc-calibration/PACS/PSF/PACSPSF_PICC-ME-TN-033_v2.2.tar.gz)
# give the most recent PACS beam profiles, and encircled energy fractions (EEFs) for different aperture sizes in the
# various scanning modes.

# PIETER:

# Zie bijlage voor de tabellen met de correcties voor de enclosed energy fraction (EEF) voor elke aperture size
# for PACS en SPIRE.
# Voor de aperture size moet ge mogelijks voor elke band een andere waarde gebruiken
# (door evt de beamsize in rekening te brengen ofzo).
# De flux in elke band moet dan gedeeld worden door de correction factor voor die band, gebruik
# makend van de aperture size in die band.

# De Growth_Curve_Final_XXmicron.dat geven de aperture size in arcsec en de correction factor voor de 2 PACS bands.

# De PSF_correction_HATLAS_SPIRE.dat geeft de aperture size in arcsec en dan de correction factors for
# 250um, 350um en 500um.

# Deze corrections zijn voor een centrale pointsource en zijn dus een soort minimum correctie voor een extended source.
# Deze minimum correcties worden doorgaands toegepast op extended sources omdat ze een goed genoege 1ste orde
# benadering zijn.

# -----------------------------------------------------------------

class PhotoMeter(PhotometryComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PhotoMeter, self).__init__(config)

        # The list of image frames
        self.images = dict()

        # The corresponding error maps
        self.errors = dict()

        # The disk ellipse
        self.disk_ellipse = None

        # The SED
        self.sed = None

        # The SEDFetcher
        self.sed_fetcher = None

        # The differences
        self.differences = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new PhotoMeter instance
        photometer = cls(arguments.config)

        # Set the modeling path
        photometer.config.path = arguments.path

        # Return the new instance
        return photometer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the truncated images
        self.load_images()

        # 3. Calculate the Enclosed Energy Fractions
        self.calculate_eefs()

        # 3. Do the photometry
        self.do_photometry()

        # 4. Get the photometric flux points from the literature for comparison
        self.get_references()

        # 5. Calculate the differences between the calculated photometry and the reference SEDs
        self.calculate_differences()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PhotoMeter, self).setup()

        # Create an observed SED
        self.sed = ObservedSED()

        # Create an SEDFetcher instance
        self.sed_fetcher = SEDFetcher()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images and error maps ...")

        # Load all data
        for name in self.dataset.names:

            # Debugging
            log.debug("Loading the data and error map for the " + name + " image ...")

            # Load the frame
            frame = self.dataset.get_frame(name)

            # Load the error map
            errors = self.dataset.get_errors(name)

            # Debugging
            log.debug("Converting the " + name + " to Jy ...")

            # CONVERSION TO JANSKY

            # Convert from MJy/sr to Jy/sr
            conversion_factor = 1.0
            conversion_factor *= 1e6

            # Conversion from Jy / sr to Jy / pixel
            pixelscale = self.images[name].average_pixelscale
            pixel_factor = (1.0 / pixelscale ** 2).to("pix2/sr").value
            conversion_factor /= pixel_factor

            # CONVERT IMAGE
            frame *= conversion_factor

            # CONVERT ERROR MAP
            errors *= conversion_factor

            # Add to the appropriate dictionary
            self.images[name] = frame
            self.errors[name] = errors

    # -----------------------------------------------------------------

    def calculate_eefs(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def do_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the photometry calculation ...")

        # The object that keeps track of PSFs and convolution kernels
        kernels = AnianoKernels()

        # Loop over all the images
        for name in self.images:

            # Debugging
            log.debug("Performing photometry for the " + name + " image ...")

            # Debugging
            log.debug("Calculating the total flux and flux error ...")

            # Calculate the total flux in Jansky
            flux = self.images[name].sum()

            # Apply correction for EEF of aperture
            if "Pacs" in name or "SPIRE" in name: flux *= self.get_aperture_correction_factor(self.images[name], kernels)

            # Calculate the total flux error in Jansky
            flux_error = self.errors[name].sum()

            # Create errorbar
            errorbar = ErrorBar(float(flux_error))

            # Add this entry to the SED
            self.sed.add_entry(self.images[name].filter, flux, errorbar)

    # -----------------------------------------------------------------

    def get_references(self):

        """
        This function ...
        :return:
        """

        # Specify which references should be consulted
        self.sed_fetcher.config.catalogs = ["GALEX", "2MASS", "SINGS", "LVL", "Spitzer", "Spitzer/IRS", "IRAS", "IRAS-FSC", "S4G", "Brown", "Planck"]

        # Fetch the reference SEDs
        self.sed_fetcher.run(self.galaxy_name)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Get list of instruments, bands and fluxes of the calculated SED
        instruments = self.sed.instruments()
        bands = self.sed.bands()
        fluxes = self.sed.fluxes(unit="Jy", add_unit=False)

        # The number of data points
        number_of_points = len(instruments)

        # Initialize data and names
        reference_labels = self.sed_fetcher.seds.keys()
        data = [[] for _ in range(len(reference_labels)+3)]
        names = ["Instrument", "Band", "Flux"]
        for label in reference_labels:
            names.append(label)

        # Loop over the different points in the calculated SED
        for i in range(number_of_points):

            # Add instrument, band and flux
            data[0].append(instruments[i])
            data[1].append(bands[i])
            data[2].append(fluxes[i])

            column_index = 3

            # Loop over the different reference SEDs
            for label in reference_labels:

                relative_difference = None

                # Loop over the data points in the reference SED
                for j in range(len(self.sed_fetcher.seds[label].table["Wavelength"])):

                    if self.sed_fetcher.seds[label].table["Instrument"][j] == instruments[i] and self.sed_fetcher.seds[label].table["Band"][j] == bands[i]:

                        difference = fluxes[i] - self.sed_fetcher.seds[label].table["Flux"][j]
                        relative_difference = difference / fluxes[i] * 100.

                        # Break because a match has been found within this reference SED
                        break

                # Add percentage to the table (or None if no match was found in this reference SED)
                data[column_index].append(relative_difference)

                column_index += 1

        # Create table of differences
        self.differences = tables.new(data, names=names)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write SED table
        self.write_sed()

        # Write the differences
        self.write_differences()

        # Plot the SED
        self.plot_sed()

        # Plot the SED with references
        self.plot_sed_with_references()

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing SED to a data file ...")

        # Save the SED
        self.sed.save(self.observed_sed_path)

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the percentual differences with reference fluxes to a data file ...")

        # Determine the full path to the output file
        path = fs.join(self.phot_path, "differences.dat")

        # Save the differences table
        tables.write(self.differences, path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def plot_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SED ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter(self.galaxy_name)

        # Add the SED
        plotter.add_observed_sed(self.sed, "PTS")

        # Determine the full path to the plot file
        path = fs.join(self.phot_path, "sed.pdf")
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_sed_with_references(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SED with reference fluxes ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter(self.galaxy_name)

        # Add the SED
        plotter.add_observed_sed(self.sed, "PTS")

        # Add the reference SEDs
        for label in self.sed_fetcher.seds: plotter.add_observed_sed(self.sed_fetcher.seds[label], label)

        # Determine the full path to the plot file
        path = fs.join(self.phot_path, "sed_with_references.pdf")
        plotter.run(path)

    # -----------------------------------------------------------------

    def get_aperture_correction_factor(self, frame, aniano):

        """
        This function ...
        :param frame:
        :param aniano:
        :return:
        """

        filter_name = str(frame.filter)

        input_dict = dict()

        # Set cutout
        input_dict["cutout"] = frame.data
        input_dict["pix_arcsec"] = frame.average_pixelscale.to("arcsec").value


        truncation_ellipse_sky = self.truncation_ellipse
        truncation_ellipse_image = truncation_ellipse_sky.to_pixel(frame.wcs)


        # PACS BLUE
        if filter_name == "Pacs blue":

            ## USE psf.data !!
            psf_path = aniano.get_psf_path(self.reference_filter) # reference filter = pacs red filter

            annulus_inner = self.sky_annulus_inner("Pacs blue")
            annulus_outer = self.sky_annulus_outer("Pacs blue")

        # PACS GREEN
        elif filter_name == "Pacs green":

            psf_path = aniano.get_psf_path(self.reference_filter) # reference filter = pacs red filter

            annulus_inner = self.sky_annulus_inner("Pacs green")
            annulus_outer = self.sky_annulus_outer("Pacs green")

        # PACS RED
        elif filter_name == "Pacs red":

            psf_path = aniano.get_psf_path(self.reference_filter)  # reference filter = pacs red filter

            annulus_inner = self.sky_annulus_inner("Pacs red")
            annulus_outer = self.sky_annulus_outer("Pacs red")

        # SPIRE PSW
        elif filter_name == "SPIRE PSW":

            psf_path = aniano.get_psf_path(frame.filter)

            annulus_inner = self.sky_annulus_inner(filter_name)
            annulus_outer = self.sky_annulus_outer(filter_name)

        # SPIRE PMW
        elif filter_name == "SPIRE PMW":

            psf_path = aniano.get_psf_path(frame.filter)

            annulus_inner = self.sky_annulus_inner(filter_name)
            annulus_outer = self.sky_annulus_outer(filter_name)

        # SPIRE PLW
        elif filter_name == "SPIRE PLW":

            psf_path = aniano.get_psf_path(frame.filter)

            annulus_inner = self.sky_annulus_inner(filter_name)
            annulus_outer = self.sky_annulus_outer(filter_name)

        # INVALID FILTER
        else: raise ValueError("Invalid filter: '" + filter_name + "'")


        # Load convolution kernel
        psf = ConvolutionKernel.from_file(psf_path)
        psf.prepare_for(frame)

        # SET INPUT DICT

        input_dict["psf"] = psf.data

        annulus_inner_factor = annulus_inner.radius.x / truncation_ellipse_image.radius.x
        annulus_outer_factor = annulus_outer.radius.x / truncation_ellipse_image.radius.x

        assert annulus_inner_factor == annulus_inner.radius.y / truncation_ellipse_image.radius.y
        assert annulus_outer_factor == annulus_outer.radius.y / truncation_ellipse_image.radius.y


        input_dict["semimaj_pix"] = truncation_ellipse_image.radius.x
        input_dict["axial_ratio"] = truncation_ellipse_image.radius.x / truncation_ellipse_image.radius.y
        input_dict["angle"] = truncation_ellipse_image.angle.to("deg").value
        input_dict["centre_i"] = truncation_ellipse_image.center.x
        input_dict["centre_j"] = truncation_ellipse_image.center.y

        input_dict["annulus_inner"] = annulus_inner_factor
        input_dict["annulus_outer"] = annulus_outer_factor


        # Calculate the aperture correction factor
        return calculate_aperture_correction(input_dict)

# -----------------------------------------------------------------

def calculate_aperture_correction(input_dict):

    """
    # Define function that uses provided beam profile to aperture-correct photometry
    Entries in input_dict:

    psf_path: Either a string giving the path to FITS file that contains the PSF, or a False boolean
              (in which case an airy disc PSF will be assumed).
    cutout: Array upon which photometry is being perfomred upon.

    pix_arcsec: The width, in arscec, of the pixels in the map photometry is being performed upon
                (this is needed in case there is a pixel size mismatch with PSF).
    semimaj_pix: Semi-major axis of photometric aperture, in pixels.
    axial_ratio: Axial ratio of photometryic aperture.
    angle: Position angle of photometric aperture, in degrees.

    centre_i: Zero-indexed, 0th-axis coordinate (equivalent to y-axis one-indexed coordinates in FITS terms)
              of centre position of photometric aperture.
    centre_j: Zero-indexed, 1st-axis coordinate (equivalent to x-axis one-indexed coordinates in FITS terms) of
              centre position of photometric aperture.

    annulus_inner: The semi-major axis of the inner edge of the background annulus, in units of the semi-major
                   axis of the source ellipse.
    annulus_outer: The semi-major axis of the outer edge of the background annulus, in units of the semi-major axis
                   of the source ellipse.
    """

    # Import things
    import ChrisFuncs
    import lmfit
    import time
    import astropy.convolution
    import scipy.ndimage
    import astropy.io.fits
    import astropy.wcs
    import astropy.modeling

    ### REBINNING, RECENTERING AND NORMALIZATION OF PSF

    """
    # If no PSF given, assume Airy disc; else extract PSF from provided file
    if (str(input_dict['psf_path'])==False) or (input_dict['psf_path'] is None):
        psf = astropy.convolution.kernels.AiryDisk2DKernel(input_dict['psf_path']).array
    else:

        # Read in PSF, and establish pixel size
        psf_in, psf_header = astropy.io.fits.getdata(input_dict['psf_path'], header=True)
        psf_wcs = astropy.wcs.WCS(psf_header)
        if psf_wcs.wcs.has_cd():
            psf_cdelt = psf_wcs.wcs.cd.max()
        else:
            psf_cdelt = psf_wcs.wcs.cdelt.max()
        psf_cdelt_arcsec = abs( psf_cdelt * 3600.0 )

        # If PSF pixel size is different to map pixel size, rescale PSF accordingly
        if (input_dict['pix_arcsec']/psf_cdelt_arcsec)>1.001 or (input_dict['pix_arcsec']/psf_cdelt_arcsec)<0.999:

            # 'Trim' the edges of the input PSF until the rescaled PSF has odd dimensions
            psf_even = True
            while psf_even:
                zoom_factor = float(psf_cdelt_arcsec) / float(input_dict['pix_arcsec'])
                psf = scipy.ndimage.zoom(psf_in, (zoom_factor,zoom_factor), mode='nearest')
                if (psf.shape[0]%2!=0) and (psf.shape[1]%2!=0):
                    psf_even = False
                else:
                    psf_in = psf_in[1:,1:]
                    psf_in = psf_in[:-1,:-1]

        # Else, if pixel sizes are already the same, leave as-is
        elif ((input_dict['pix_arcsec']/psf_cdelt_arcsec)>=0.999) and ((input_dict['pix_arcsec']/psf_cdelt_arcsec)<=0.001):
            psf = psf_in.copy()

    # Normalise PSF
    psf /= np.nansum(psf)

    """

    psf = input_dict["psf"]

    #####

    # Extract cutout
    cutout = input_dict["cutout"]

    # Produce mask for pixels we care about for fitting (ie, are inside photometric aperture and background annulus)
    mask = ChrisFuncs.Photom.EllipseMask(cutout, input_dict['semimaj_pix'], input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j']) #*band_dict['annulus_outer']

    ##
    from pts.magic.tools import plotting
    plotting.plot_mask(mask)
    ##

    # Produce guess values
    initial_sersic_amplitide = cutout[ input_dict['centre_i'], input_dict['centre_j'] ]
    initial_sersic_r_eff = input_dict['semimaj_pix'] / 10.0
    initial_sersic_n = 1.0
    initial_sersic_x_0 = input_dict['centre_j']
    initial_sersic_y_0 = input_dict['centre_i']
    initial_sersic_ellip = ( input_dict['axial_ratio'] - 1.0 ) / input_dict['axial_ratio']
    initial_sersic_theta = np.deg2rad( input_dict['angle'] )

    # Produce sersic model from guess parameters, for time trials
    sersic_x, sersic_y = np.meshgrid( np.arange(cutout.shape[1]), np.arange(cutout.shape[0]) )
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=initial_sersic_amplitide, r_eff=initial_sersic_r_eff, n=initial_sersic_n, x_0=initial_sersic_x_0, y_0=initial_sersic_y_0, ellip=initial_sersic_ellip, theta=initial_sersic_theta)
    sersic_map = sersic_model(sersic_x,sersic_y)

    # Make sure that PSF array is smaller than sersic model array (as required for convolution); if not, remove its edges such that it is
    if psf.shape[0] > sersic_map.shape[0] or psf.shape[1] > sersic_map.shape[1]:

        excess = max( psf.shape[0]-sersic_map.shape[0], psf.shape[1] - sersic_map.shape[1] )
        border = int( np.round( np.ceil( float(excess) / 2.0 ) - 1.0 ) )
        psf = psf[border:,border:]
        psf = psf[:-border,:-border]


    # Determine whether FFT convolution or direct convolution is faster for this kernel,
    # using sersic model produced with guess parameters
    time_fft = time.time()
    conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    time_fft = time.time() - time_fft
    time_direct = time.time()
    conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)
    time_direct = time.time() - time_direct

    if time_fft < time_direct: use_fft = True
    else: use_fft = False


    # Set up parameters to fit galaxy with 2-dimensional sersic profile
    params = lmfit.Parameters()
    params.add('sersic_amplitide', value=initial_sersic_amplitide, vary=True)
    params.add('sersic_r_eff', value=initial_sersic_r_eff, vary=True, min=0.0, max=input_dict['semimaj_pix'])
    params.add('sersic_n', value=initial_sersic_n, vary=True, min=0.1, max=10)
    params.add('sersic_x_0', value=initial_sersic_x_0, vary=False)
    params.add('sersic_y_0', value=initial_sersic_y_0, vary=False)
    params.add('sersic_ellip', value=initial_sersic_ellip, vary=True, min=0.5*initial_sersic_ellip, max=0.5*(1.0-initial_sersic_ellip)+initial_sersic_ellip)
    params.add('sersic_theta', value=initial_sersic_theta, vary=False)

    # Solve with LMfit to find parameters of best-fit sersic profile
    result = lmfit.minimize(chi_squared_sersic, params, args=(cutout, psf, mask, use_fft), method='leastsq', ftol=1E-5, xtol=1E-5)

    # Extract best-fit results
    sersic_amplitide = result.params['sersic_amplitide'].value
    sersic_r_eff = result.params['sersic_r_eff'].value
    sersic_n = result.params['sersic_n'].value
    sersic_x_0 = result.params['sersic_x_0'].value
    sersic_y_0 = result.params['sersic_y_0'].value
    sersic_ellip = result.params['sersic_ellip'].value
    sersic_theta = result.params['sersic_theta'].value

    # Construct underlying sersic map and convolved sersic map, using best-fit parameters
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
    sersic_map = sersic_model(sersic_x, sersic_y)

    if use_fft: conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    else: conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)



    # Determine annulus properties before proceeding with photometry
    bg_inner_semimaj_pix = input_dict['semimaj_pix'] * input_dict['annulus_inner']
    bg_width = (input_dict['semimaj_pix'] * input_dict['annulus_outer']) - bg_inner_semimaj_pix

    # Evaluate pixels in source aperture and background annulus  unconvoled sersic map
    sersic_ap_calc = ChrisFuncs.Photom.EllipseSum(sersic_map, input_dict['semimaj_pix'], input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])
    sersic_bg_calc = ChrisFuncs.Photom.AnnulusSum(sersic_map, bg_inner_semimaj_pix, bg_width, input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])

    # Background-subtract and measure unconvoled sersic source flux
    sersic_bg_clip = ChrisFuncs.SigmaClip(sersic_bg_calc[2], median=False, sigma_thresh=3.0)
    sersic_bg_avg = sersic_bg_clip[1]
    sersic_ap_sum = sersic_ap_calc[0] - (sersic_ap_calc[1] * sersic_bg_avg)

    # Evaluate pixels in source aperture and background annulus in convolved sersic map
    conv_ap_calc = ChrisFuncs.Photom.EllipseSum(conv_map, input_dict['semimaj_pix'], input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])
    conv_bg_calc = ChrisFuncs.Photom.AnnulusSum(conv_map, bg_inner_semimaj_pix, bg_width, input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])

    # Background-subtract and measure convolved sersic source flux
    conv_bg_clip = ChrisFuncs.SigmaClip(conv_bg_calc[2], median=False, sigma_thresh=3.0)
    conv_bg_avg = conv_bg_clip[1]
    conv_ap_sum = conv_ap_calc[0] - (conv_ap_calc[1] * conv_bg_avg)



    # Find difference between flux measued on convoled and unconvoled sersic maps
    ap_correction = np.nanmax([ 1.0, (sersic_ap_sum/conv_ap_sum) ])

    # Return aperture correction
    return ap_correction

# -----------------------------------------------------------------

def chi_squared_sersic(params, cutout, psf, mask, use_fft, lmfit=True):

    """
    This function defines LMfit convolved-sersic function
    :param params:
    :param cutout:
    :param psf:
    :param mask:
    :param use_fft:
    :param lmfit:
    :return:
    """

    # Import things
    import astropy.modeling
    import astropy.convolution


    # Extract variable parameters
    sersic_amplitide = params['sersic_amplitide'].value
    sersic_r_eff = params['sersic_r_eff'].value
    sersic_n = params['sersic_n'].value
    sersic_x_0 = params['sersic_x_0'].value
    sersic_y_0 = params['sersic_y_0'].value
    sersic_ellip = params['sersic_ellip'].value
    sersic_theta = params['sersic_theta'].value

    # Generate sersic model given current parameters
    sersic_x, sersic_y = np.meshgrid( np.arange(cutout.shape[1]), np.arange(cutout.shape[0]) )
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
    sersic_map = sersic_model( sersic_x, sersic_y )

    # Convolve sersic model with PSF
    if use_fft: conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    else: conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)

    # Calculate residuals, filtered by mask
    residuals = cutout - conv_map
    mask[ np.where(mask==0.0) ] = np.nan
    residuals *= mask
    residuals.flatten()
    residuals = residuals[ np.where( np.isnan(residuals)==False ) ]

    # Return residuals
    if lmfit: return residuals**2.0
    else: return residuals, cutout-conv_map

# -----------------------------------------------------------------

#input_dict = {}
#input_dict['psf_path'] = '/home/herdata/spx7cjc/Beams/SPIRE_250.fits'
##input_dict['cutout'] = astropy.io.fits.getdata('/home/saruman/spx7cjc/DustPedia/SPIRE/Cutouts/DustPedia/NGC4030_SPIRE_250.fits')
#input_dict["cutout"] = astropy.io.fits.getdata("NGC4030_SPIRE_250.fits")
#input_dict['pix_arcsec'] = 6.0
#input_dict['semimaj_pix'] = 41.0
#input_dict['axial_ratio'] = 1.1795263352195566
#input_dict['angle'] = 115.16660752050387
#input_dict['centre_i'] = 300.0
#input_dict['centre_j'] = 300.0
#input_dict['annulus_inner'] = 1.25
#input_dict['annulus_outer'] = 1.601

#ap_corr = StandaloneApCorrect(input_dict)

# -----------------------------------------------------------------
