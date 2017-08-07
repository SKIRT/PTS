#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.makeinfofile Creating information files with statistics on EAGLE SKIRT-runs.
#
# The main function in this module creates an information file with statistics on the results of an EAGLE SKIRT-run.
# The information file has a simple text format and includes statistics on the data exported from the EAGLE snapshot
# (i.e. the input data for the SKIRT simulation), on the SKIRT setup (such as the dust grid), and on the results
# of the SKIRT simulation (such as fluxes in various bands).

# -----------------------------------------------------------------

# Import standard modules
import os.path
import numpy as np
import pyfits
from scipy.ndimage.filters import gaussian_filter

# Import the relevant PTS classes and modules
from ..core.tools import archive as arch
from ..core.filter.broad import BroadBandFilter

# -----------------------------------------------------------------

# private list of standard filters for which integrated fluxes should be calculated
_filterspecs = ( \
                "GALEX.FUV","GALEX.NUV",
                "SDSS.u","SDSS.g","SDSS.r","SDSS.i","SDSS.z",
                "UKIDSS.Z","UKIDSS.Y","UKIDSS.J","UKIDSS.H","UKIDSS.K",
                "Johnson.U","Johnson.B","Johnson.V","Johnson.R","Johnson.I","Johnson.J","Johnson.M",
                "2MASS.J", "2MASS.H", "2MASS.K",
                "IRAS.12","IRAS.25","IRAS.60","IRAS.100",
                "IRAC.I1","IRAC.I2","IRAC.I3","IRAC.I4",
                "MIPS.24", "MIPS.70", "MIPS.160",
                "WISE.W1", "WISE.W2", "WISE.W3", "WISE.W4",
                "Pacs.blue","Pacs.green","Pacs.red",
                "SPIRE.PSW","SPIRE.PMW","SPIRE.PLW",
                "SCUBA2.450", "SCUBA2.850"
                )

# private list of uniform filters for which integrated fluxes should be calculated;
# for each filter, specify the filter spec and the wavelength range in micron)
_uniformfilters = ( ("Uniform_8_1000",8,1000), ("Uniform_3_1100",3,1100))

# private global dictionary to hold filter objects: (key,value) = (filterspec, Filter object)
_filters = { }

## This private function ensures that the filters relevant for this module are loaded into the private global
# dictionary \em _filters, using the filter specification string as a key and the corresponding Filter instance
# as a value. When the function is called for the first time, it creates the Filter objects and stores them
# in the dictionary. Once the dictionary has been loaded, subsequent invocations of the function don't do anything.
def _loadfilters():
    if len(_filters) == 0:
        # standard filters
        for filterspec in _filterspecs:
            _filters[filterspec] = BroadBandFilter(filterspec)
        # uniform filters
        for filterspec,wavemin,wavemax in _uniformfilters:
            _filters[filterspec] = BroadBandFilter((wavemin,wavemax))

# -----------------------------------------------------------------

## This function creates an information file with statistics on the results of an EAGLE SKIRT-run.
# The information file has a simple text format and includes statistics on the data exported from the EAGLE snapshot
# (i.e. the input data for the SKIRT simulation), on the SKIRT setup (such as the dust grid), and on the results
# of the SKIRT simulation (such as fluxes in various bands).
#
# The information file is placed in the SKIRT-run visualization directory, and is named "prefix_info.txt".
#
def makeinfofile(skirtrun, snaptag):
    simulation = skirtrun.simulation()

    # get the redshift corresponding to the snapshot tag
    redshift = { 28:   0.0000,
                 27:   0.1006,
                 26:   0.1827,
                 25:   0.2709,
                 24:   0.3657,
                 23:   0.5031,
                 22:   0.6152,
                 21:   0.7356,
                 20:   0.8651,
                 19:   1.0041,
                 18:   1.2593,
                 17:   1.4867,
                 16:   1.7370,
                 15:   2.0124,
                 14:   2.2370,
                 13:   2.4784,
                 12:   3.0165,
                 11:   3.5280,
                 10:   3.9837,
                  9:   4.4852,
                  8:   5.0372,
                  7:   5.4874,
                  6:   5.9712,
                  5:   7.0496,
                  4:   8.0746,
                  3:   8.9879,
                  2:   9.9930,
                  1:  15.1323,
                  0:  20.0000 } [snaptag]

    # create the info dict that will be saved at the end
    info = { }
    info["skirt_run_id"] = skirtrun.runid()
    info["original_redshift"] = redshift

    # load statistics on the EAGLE data from the info file in the input folder
    inpath = skirtrun.inpath()
    infile = arch.listdir(inpath, "_info.txt")[0]
    for line in arch.opentext(os.path.join(inpath,infile)):
        if not line.startswith("#"):
            segments = line.split()
            info[segments[0]] = float(segments[2])

    # compute the dust mass from the gas input file and the parameters in the ski file
    infile = arch.listdir(inpath, "_gas.dat")[0]
    info["exported_mass_cold_gas"], info["exported_mass_negative_cold_gas"], \
    info["exported_mass_metallic_gas"], info["exported_mass_negative_metallic_gas"] = \
            loadgasmasses(os.path.join(inpath,infile), simulation.maximumtemperature())
    info["exported_mass_dust"] = info["exported_mass_metallic_gas"] * simulation.dustfraction()
    info["exported_mass_negative_dust"] = info["exported_mass_negative_metallic_gas"] * simulation.dustfraction()

    # gather SKIRT setup statistics
    info["setup_mass_dust"] = simulation.dustmass()
    info["setup_mass_dust_grid"] = simulation.dustgridmass()
    info["setup_mass_cold_gas"] = simulation.coldgasmass()
    info["setup_particles_cold_gas"] = simulation.coldgasparticles()
    info["setup_mass_metallic_gas"] = simulation.metallicgasmass()
    info["setup_initial_mass_stars"] = simulation.initialstellarmass()
    info["setup_mass_hii_regions"] = simulation.hiiregionmass()
    info["setup_luminosity_stars"] = simulation.stellarluminosity()
    info["setup_luminosity_hii_regions"] = simulation.hiiregionluminosity()
    info["setup_cells_dust_grid"], info["setup_optical_depth_maximum"], \
        info["setup_optical_depth_percentile90"] = map(float,simulation.dustcellstats())
    info["setup_distance_instrument"] = simulation.instrumentdistance()

    # gather SKIRT operating statistics
    info["skirt_run_wall_time"], info["skirt_run_peak_memory"], \
    info["skirt_run_self_absorption_cycles"] = operatingStatistics(simulation)

    # load filters and wavelength grid
    _loadfilters()
    wavelengths = simulation.wavelengths()
    shifted_wavelengths = wavelengths * (1 + redshift)

    # create a mask that removes the carbon line emission peaks from the dust continuum emission
    cmask = (np.abs(wavelengths-157.5)>3) & (np.abs(wavelengths-360)>20) & (np.abs(wavelengths-600)>20)

    # define properties for the Herschel instruments used in determining dust temperature
    h_filterspecs = ("Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW")
    # --> beam FWHM & area for PACS from Herschel PACS observer's manual, July 2013
    # --> beam FWHM & area for SPIRE from Ciesla et al. 2012, A&A 543, A161
    h_beam_fwhms = (12, 18.2, 24.5, 36.0)       # in arcsec
    h_beam_areas = (180, 423, 751, 1587)        # in arcsec^2
    # --> flux limit for PACS = noise level from Cortese et al. 2014
    # --> flux limit for SPIRE = confusion noise from Nguyen et al. 2010 A&A 518, L5
    h_flux_limits = (5, 5.8, 6.3, 6.8)          # in mJy/beam
    # --> convert to MJy/sr
    h_flux_limits = np.array(h_flux_limits) / np.array(h_beam_areas) * (648000/np.pi)**2 * 1e-9

    # wavelength ranges and masks for the H-alpha emission line and the continuum immediately next to it
    wline0 = 0.6565
    wline1, wline2 = 0.6520, 0.6610
    wcont1, wcont2 = 0.6420, 0.6700   # upper limit should be 0.6650 for more accurate results with finer grid
    linemask = ((wavelengths>=wline1) & (wavelengths<=wline2))
    contmask = ((wavelengths>=wcont1) & (wavelengths<=wline1)) | ((wavelengths>=wline2) & (wavelengths<=wcont2))

    # create uniform filter over the wavelength range of the H-alpha emission peak
    halpha_filter = BroadBandFilter((0.6500,0.6631))     # limits so that pivot wavelength coincides with peak
    # create a mask that removes the H-alpha emission peak from the continuum emission
    halpha_mask = np.abs(wavelengths-0.6565)>0.0010

    # gather statistics on fluxes received by each instrument
    for name in simulation.instrumentnames():
        # maximum flux in Jy
        fluxdensities = simulation.fluxdensities(name, unit='Jy')
        info["instr_"+name+"_fluxdensity_maximum"] = fluxdensities.max()

        # for each filter, calculate integrated flux density and absolute magnitude
        for filterspec,filterobject in _filters.iteritems():
            fluxdensity = regularfluxdensity(simulation, name, [1], wavelengths, None, filterobject)
            addfluxinfo(info, simulation, name, filterspec, fluxdensity)
            if redshift>0.01:
                fluxdensity = regularfluxdensity(simulation, name, [1], shifted_wavelengths, None, filterobject)
                fluxdensity = distantflux(fluxdensity, simulation.instrumentdistance(unit='m'), redshift)
            addfluxinfo(info, simulation, name, filterspec, fluxdensity, "observer")

        # for the Herschel filters, calculate flux and magnitude excluding the carbon line emission peaks
        for filterspec in ("Pacs.blue","Pacs.green","Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW","SCUBA2.450","SCUBA2.850"):
            fluxdensity = regularfluxdensity(simulation, name, [1], wavelengths, cmask, _filters[filterspec])
            addfluxinfo(info, simulation, name, filterspec, fluxdensity, "continuum")
            if redshift>0.01:
                fluxdensity = regularfluxdensity(simulation, name, [1], shifted_wavelengths, cmask, _filters[filterspec])
                fluxdensity = distantflux(fluxdensity, simulation.instrumentdistance(unit='m'), redshift)
            addfluxinfo(info, simulation, name, filterspec, fluxdensity, "continuum_observer")
            fluxdensity = regularfluxdensity(simulation, name, [2,3], wavelengths, cmask, _filters[filterspec])
            addfluxinfo(info, simulation, name, filterspec, fluxdensity, "hii_continuum")
            fluxdensity = regularfluxdensity(simulation, name, [4], wavelengths, cmask, _filters[filterspec])
            addfluxinfo(info, simulation, name, filterspec, fluxdensity, "other_continuum")

        # for the Herschel filters used in determining dust temperature, calculate the "limited" flux and magnitude
        # (i.e. ignoring pixels with a flux under a specific limit, and still excluding the carbon line emission peaks)
        if (simulation.nsimpleinstruments() + simulation.nfullinstruments()) > 0:
            for filterspec,fwhm,fluxlimit in zip(h_filterspecs, h_beam_fwhms, h_flux_limits):
                fluxdensity = limitedfluxdensity(simulation, name, wavelengths, cmask, _filters[filterspec], fwhm, fluxlimit)
                addfluxinfo(info, simulation, name, filterspec, fluxdensity, "limited")

        # calculate the H-alpha flux density
        fluxdensities = simulation.fluxdensities(name, unit='W/m2/micron')
        rico, f0 = np.polyfit(wavelengths[contmask], fluxdensities[contmask], 1)
        linefluxdensities = fluxdensities - (f0 + rico*wavelengths)
        linefluxdensities[linefluxdensities<0] = 0
        linefluxdensity = np.trapz(x=wavelengths[linemask], y=linefluxdensities[linemask]) / wline0
        linefluxdensity = simulation.convert(linefluxdensity, from_unit='W/m2/micron', to_unit='Jy', wavelength=wline0)
        addfluxinfo(info, simulation, name, "halpha", linefluxdensity)

    # save the info file
    infofilepath = os.path.join(skirtrun.vispath(), simulation.prefix() + "_info.txt")
    infofile = open(infofilepath, 'w')
    infofile.write('# Information file for SKIRT-run {}\n'.format(skirtrun.runid()))
    infofile.write('# cells : 1\n')
    infofile.write('# particles : 1\n')
    infofile.write('# cycles : 1\n')
    infofile.write('# luminosity : Lsun\n')
    infofile.write('# mass : Msun\n')
    infofile.write('# distance : pc\n')
    infofile.write('# fluxdensity : Jy\n')
    infofile.write('# magnitude : mag\n')
    infofile.write('# memory : GB\n')
    infofile.write('# time : s\n')
    maxkeylen = max(map(len,info.keys()))
    for key in sorted(info.keys()):
        valueformat = ".0f" if "_particles_" in key or "_cells_" in key or "_id" in key or "_cycles" in key else ".9e"
        infofile.write( ("{0:"+str(maxkeylen)+"} = {1:16"+valueformat+"}\n").format(key, info[key]) )
    infofile.close()

    # report success
    print "Created info file " + infofilepath

# -----------------------------------------------------------------

## This function calculates and returns the positive and negative cold and metallic gas masses
# from the specified gas input file and maximum temperature; it returns a tuple
# (positive cold gas mass, negative cold gas mass, positive metallic gas mass, negative metallic gas mass)

def loadgasmasses(inpath, Tmax):
    data = np.loadtxt(arch.opentext(inpath), usecols=(4,5,6), unpack=True)
    if len(data)>0:
        M,Z,T = data
        mask_pos = (T<=Tmax) & (M>0)
        mask_neg = (T<=Tmax) & (M<0)
        return M[mask_pos].sum(), -M[mask_neg].sum(), (M[mask_pos]*Z[mask_pos]).sum(), -(M[mask_neg]*Z[mask_neg]).sum()
    else:
        return 0,0,0,0

# -----------------------------------------------------------------

## This function adds the specified flux density and the corresponding magnitude to the info dictionary.
# The function takes the following arguments:
#  - \em info: the info dictionary in which to store the information.
#  - \em simulation: a SkirtSimulation instance representing the simulation for which to perform the action.
#  - \em instrumentname: the name of the instrument; for use in constructing the info dictionary key.
#  - \em filterspec: the string that specifies the filter; for use in constructing the info dictionary key.
#  - \em fluxdensity: the flux density to be added to the dictionary, in Jy. If None, the function does nothing.
#  - \em fluxtype: a string specifying the flux type (no underscore); for use in constructing the info dictionary key.
#
def addfluxinfo(info, simulation, instrumentname, filterspec, fluxdensity, fluxtype=""):
    if fluxdensity is not None:
        magnitude = simulation.absolutemagnitude(fluxdensity, simulation.instrumentdistance(unit='pc'),
                                                 fluxdensity_unit='Jy', distance_unit='pc')
        filtername = filterspec.replace(".","_").lower()
        if fluxtype!="": fluxtype = "_"+fluxtype
        info["instr_"+instrumentname+"_fluxdensity_"+filtername+fluxtype] = fluxdensity
        info["instr_"+instrumentname+"_magnitude_"+filtername+fluxtype] = magnitude

# -----------------------------------------------------------------

## This function calculates and returns the flux density (in Jy) in a given band based on the total sed or on some of
# the component seds of a particular instrument. The function convolves the sed over the wavelengths of the specified
# filter to obtain an average value.
#
# The function takes the following arguments:
#  - \em simulation: a SkirtSimulation instance representing the simulation for which to perform the calculation.
#  - \em instrumentname: the name of the instrument (as listed in the ski file) for which to perform the calculation.
#  - \em columns: a sequence of zero-based column indices in the sed data file; these flux densities are summed.
#  - \em wavelengths: the simulation's wavelength grid (this could be retrieved from the simulation but happens to
#                     be available at the caller site already).
#  - \em cmask: a mask indicating the wavelengths to include/exclude in the convolution with the filter, or None.
#  - \em filterobject: the filter defining the band for which to perform the convolution.
#
def regularfluxdensity(simulation, instrumentname, columns, wavelengths, cmask, filterobject):

    # get the path for the sed data file corresponding to this instrument
    sedpaths = filter(lambda fn: ("_"+instrumentname+"_") in fn, simulation.seddatpaths())
    if len(sedpaths) != 1: return None

    # we do this inside a try block in case the specified columns are not available
    try:
        # load the flux densities and convert them to a per-wavelength unit, as required by Filter.convolve()
        fluxdensities = np.loadtxt(arch.opentext(sedpaths[0]), usecols=columns, unpack=True)
        if len(columns)>1: fluxdensities = fluxdensities.sum(axis=0)
        fluxdensities = simulation.convert(fluxdensities, to_unit='W/m2/micron', quantity='fluxdensity',
                                           wavelength=wavelengths)

        # perform the convolution, using the mask if requested
        if cmask is None:
            fluxdensity = filterobject.convolve(wavelengths, fluxdensities)
        else:
            fluxdensity = filterobject.convolve(wavelengths[cmask], fluxdensities[cmask])

        # convert the final flux density to Jy using the filter's pivot wavelength
        fluxdensity = simulation.convert(fluxdensity, from_unit='W/m2/micron', to_unit='Jy',
                                         wavelength=filterobject.pivotwavelength())
        return fluxdensity
    except:
        return None

# -----------------------------------------------------------------

## This function calculates and returns the "limited" flux density (in Jy) in a given band
# based on the data cube of a particular instrument. The function first convolves the data cube over
# the wavelengths of the specified filter to obtain an averaged image, then convolves the image with
# a Gaussion PSF, and finally integrates over the image ignoring any pixels with a value under the
# specified surface brightness limit.
#
# The function takes the following arguments:
#  - \em simulation: a SkirtSimulation instance representing the simulation for which to perform the calculation.
#  - \em instrumentname: the name of the instrument (as listed in the ski file) for which to perform the calculation.
#  - \em wavelengths: the simulation's wavelength grid (this could be retrieved from the simulation but happens to
#                     be available at the caller site already).
#  - \em cmask: a mask indicating the wavelengths to include/exclude in the convolution with the filter.
#  - \em filterobject: the filter defining the band for which to perform the convolution.
#  - \em fwhm: the FWHM of the Gaussion filter to be applied to the image, in arcsecs.
#  - \em fluxlimit: the min. surface brightness for a pixel to be included in the integration over the image, in MJy/sr.
#
def limitedfluxdensity(simulation, instrumentname, wavelengths, cmask, filterobject, fwhm, fluxlimit):

    # get the path for the data cube corresponding to this instrument
    fitspaths = filter(lambda fn: ("_"+instrumentname+"_") in fn, simulation.totalfitspaths())
    if len(fitspaths) != 1: return None

    # get the data cube and convert it to per-wavelength units
    cube = pyfits.getdata(arch.openbinary(fitspaths[0])).T
    cube = simulation.convert(cube, to_unit='W/m3/sr', quantity='surfacebrightness', wavelength=wavelengths)

    # convolve the data cube to a single frame, and convert back to per-frequency units
    frame = filterobject.convolve(wavelengths[cmask], cube[:,:,cmask])
    frame = simulation.convert(frame, from_unit='W/m3/sr', to_unit='MJy/sr', wavelength=filterobject.pivotwavelength())

    # get information on the simulated pixels (assume all frames are identical and square)
    sim_pixels = simulation.instrumentshape()[0]
    sim_pixelarea = simulation.angularpixelarea()   # in sr
    sim_pixelwidth = np.sqrt(sim_pixelarea) * 648000 / np.pi   # in arcsec

    # convolve the frame with a Gaussian of the appropriate size
    frame = gaussian_filter(frame, sigma=fwhm/sim_pixelwidth/2.35482, mode='constant')

    # get information on the observational instrument's pixels (assume pixel width ~ fwhm/3)
    obs_pixelwidth = fwhm/3   # in arcsec

    # calculate the bin size to obtain simulated pixels similar to observed pixels
    bin_pixels = find_nearest_divisor(sim_pixels, obs_pixelwidth/sim_pixelwidth)
    bin_pixelarea = sim_pixelarea * bin_pixels**2

    # rebin the frame
    frame = frame.reshape((sim_pixels//bin_pixels,bin_pixels,sim_pixels//bin_pixels,bin_pixels)).mean(axis=3).mean(axis=1)

    # integrate over the frame to obtain the total flux density and convert from MJy to Jy
    fluxdensity = frame[frame>fluxlimit].sum() * bin_pixelarea * 1e6
    return fluxdensity

## This function returns the integer divisor of the first (integer) argument that is nearest to the second (float) argument
def find_nearest_divisor(n, d):
    divisors = np.asarray(sorted(set( x for tup in ((i, n//i) for i in range(1, int(n**0.5)+1) if n % i == 0) for x in tup )))
    index = np.argmin(np.abs(divisors-d))
    return divisors[index]

# -----------------------------------------------------------------

## This helper function is part of the calculation of the luminosity distance for given redshift using the
# approximation presented by Adachi & Kasai 2012.
def T(x):
    b1 = 2.64086441
    b2 = 0.883044401
    b3 = 0.0531249537
    c1 = 1.39186078
    c2 = 0.512094674
    c3 = 0.0394382061
    x3 = x*x*x
    return np.sqrt(x) * (2.0+b1*x3+b2*x3*x3+b3*x3*x3*x3) / (1.0+c1*x3+c2*x3*x3+c3*x3*x3*x3)

## This function returns the luminosity distance in m for a given redshift using the
# approximation presented by Adachi & Kasai 2012.
def DL(z):
    c = 299792458.                  # in m/s
    Mpc = 3.08567758e22             # in m
    km = 1e3                        # in m
    H0 = 0.6777 * 100*km/Mpc        # in 1/s
    Omegam = 0.307
    s = ((1.0-Omegam)/Omegam) ** (1.0/3.0)
    return c*(1+z)/H0/((s*Omegam)**0.5)*(T(s)-T(s/(1+z)))

## This function corrects the given redshifted flux f in Jy, calculated for a local observation at given distance in m,
# for the luminosity distance corresponding to the given redshift.
def distantflux(f, d, z):
    # expect flux in Jy, instrument distance in m, redshift ; return corrected flux in Jy
    # see e.g. https://ned.ipac.caltech.edu/level5/Hogg/Hogg7.html
    return f * (d/DL(z))**2 * (1+z)

# -----------------------------------------------------------------

## This function retrieves some operating statistics from the simulation's log file. It returns a tuple as follows:
# (simulation walltime in seconds, peak memory usage in GB, number of dust self-absorption cycles).
def operatingStatistics(simulation):
    lines = open(simulation.logfilepath()).readlines()

    walltime = 0
    for line in reversed(lines):
        if " Finished simulation " in line:
            segments = line.split()
            timeindex = segments.index("s") if "s" in segments else segments.index("s.")
            walltime = float(segments[timeindex-1])
            break

    memory = 0
    for line in reversed(lines):
        if " Available memory: " in line:
            segments = line.split()
            memory = float(segments[segments.index("usage:")+1])
            break

    cycles = 0
    for line in reversed(lines):
        if " Finished " in line and " self-absorption cycle " in line:
            segments = line.split()
            cycles = int(segments[segments.index("cycle")+1])
            break

    return (walltime, memory, cycles)

# -----------------------------------------------------------------
