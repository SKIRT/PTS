#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.celldustheating Contains the CellDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import subprocess
import rpy2
import rpy2.robjects as ro
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ....core.tools import tables, inspection

# -----------------------------------------------------------------

class CellDustHeatingAnalyser(DustHeatingAnalysisComponent):
    
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
        super(CellDustHeatingAnalyser, self).__init__(config)

        # The wavelength grid used for the simulations
        self.wavelength_grid = None

        # The number of wavelengths
        self.number_of_wavelengths = None

        # The table with the absorbed luminosities
        self.absorption_table = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new HeatingAnalyser instance
        analyser = cls()

        # Set the modeling path
        analyser.config.path = arguments.path

        # Return the new instance
        return analyser

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Load the wavelength grid
        self.load_wavelength_grid()

        # Load the cell properties
        self.load_cell_properties()

        # Load the ISRF data
        #self.load_isrf()

        # Load the absorption data
        self.load_absorption()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(CellDustHeatingAnalyser, self).setup()

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grid file produced by SKIRT ...")

        # Determine the path to the wavelength grid file in heating/total
        wavelengths_path = fs.join(self.output_paths["total"], self.galaxy_name + "_wavelengths.dat")

        # Load the wavelength grid as a table
        self.wavelength_grid = tables.from_file(wavelengths_path, format="ascii")

        # Set the column names and units
        self.wavelength_grid.rename_column("col1", "Wavelength")
        self.wavelength_grid.rename_column("col2", "Delta")
        self.wavelength_grid["Wavelength"].unit = "micron"
        self.wavelength_grid["Delta"].unit = "micron"

        # Determine the number of wavelengths
        self.number_of_wavelengths = len(self.wavelength_grid)

    # -----------------------------------------------------------------

    def load_cell_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the suer
        log.info("Loading the dust cell properties ...")

        # M81_ds_cellprops.dat

        # column 1: volume (pc3)
        # column 2: density (Msun/pc3)
        # column 3: mass fraction
        # column 4: optical depth

        # Determine the path to the cell properties file in heating/total
        properties_path = fs.join(self.output_paths["total"], self.galaxy_name + "_ds_cellprops.dat")

        cell, x, y, z, J = np.loadtxt(isrf_path, usecols=column_indices, unpack=True)

    # -----------------------------------------------------------------

    def load_isrf(self):

        """
        This function ...
        :return:
        """

        # M81_ds_isrf.dat

        # Mean field intensities for all dust cells with nonzero absorption
        # column 1: dust cell index
        # column 2: x coordinate of cell center (pc)
        # column 3: y coordinate of cell center (pc)
        # column 4: z coordinate of cell center (pc)
        # column 5: J_lambda (W/m3/sr) for lambda = 0.1 micron
        # column 6: J_lambda (W/m3/sr) for lambda = 0.121153 micron
        # column 7: J_lambda (W/m3/sr) for lambda = 0.14678 micron
        # column 8: J_lambda (W/m3/sr) for lambda = 0.177828 micron
        # column 9: J_lambda (W/m3/sr) for lambda = 0.215443 micron
        # column 10: J_lambda (W/m3/sr) for lambda = 0.261016 micron
        # column 11: J_lambda (W/m3/sr) for lambda = 0.316228 micron
        # column 12: J_lambda (W/m3/sr) for lambda = 0.383119 micron
        # column 13: J_lambda (W/m3/sr) for lambda = 0.464159 micron
        # column 14: J_lambda (W/m3/sr) for lambda = 0.562341 micron
        # column 15: J_lambda (W/m3/sr) for lambda = 0.681292 micron
        # column 16: J_lambda (W/m3/sr) for lambda = 0.825404 micron
        # column 17: J_lambda (W/m3/sr) for lambda = 1 micron
        # column 18: J_lambda (W/m3/sr) for lambda = 1.21153 micron
        # column 19: J_lambda (W/m3/sr) for lambda = 1.4678 micron
        # column 20: J_lambda (W/m3/sr) for lambda = 1.77828 micron
        # column 21: J_lambda (W/m3/sr) for lambda = 2.15443 micron
        # column 22: J_lambda (W/m3/sr) for lambda = 2.61016 micron
        # column 23: J_lambda (W/m3/sr) for lambda = 3.16228 micron
        # column 24: J_lambda (W/m3/sr) for lambda = 3.83119 micron
        # column 25: J_lambda (W/m3/sr) for lambda = 4.64159 micron
        # column 26: J_lambda (W/m3/sr) for lambda = 5.62341 micron
        # column 27: J_lambda (W/m3/sr) for lambda = 6.81292 micron
        # column 28: J_lambda (W/m3/sr) for lambda = 8.25404 micron
        # column 29: J_lambda (W/m3/sr) for lambda = 10 micron

        # Determine the indices of the columns that have to be imported
        column_indices = [0,1,2,3]
        for i in range(4, 4 + self.number_of_wavelengths): column_indices.append(i)

        # Loop over the different contributions
        for contribution in self.contributions:

            # Determine the path to the output directory of the simulation
            output_path = self.output_paths[contribution]

            # Determine the path to the ISFR data file
            isrf_path = fs.join(output_path, self.galaxy_name + "_ds_isrf.dat")

            # Load the ISRF file
            columns = np.loadtxt(isrf_path, usecols=column_indices, unpack=True)

            # columns[0]: dust cell index
            # columns[1]: x coordinate of cell center
            # columns[2]: y coordinate of cell center
            # columns[3]: z coordinate of cell center
            # columns[4 -> 4 + (nwavelengths - 1)]: J_lambda

            # Integrate over the J_lambda values to get the total bolometric absorbed luminosity per cell
            luminosities = self.integrate_over_wavelengths(columns[4:4+self.number_of_wavelengths])

            #IDtot, x, y, z, Ltot = np.loadtxt(totISRFfile, usecols=(0, 1, 2, 3, 4,), unpack=True)
            #IDold, Lold = np.loadtxt(oldISRFfile, usecols=(0, 4,), unpack=True)
            #IDyng, Lyng = np.loadtxt(yngISRFfile, usecols=(0, 4,), unpack=True)
            #IDnew, Lnew = np.loadtxt(newISRFfile, usecols=(0, 4,), unpack=True)

            # write out
            #writeLabsTot('Labs_tot.dat', IDtot, x, y, z, Ltot)
            #writeLabs('Labs_old.dat', IDold, Lold)
            #writeLabs('Labs_yng.dat', IDyng, Lyng)
            #writeLabs('Labs_new.dat', IDnew, Lnew)

        def writeLabsTot(file, ID, x, y, z, Labs):

            with open(file, 'w') as f:
                f.write('# ID    x(pc)    y(pc)    z(pc)    Labs(W) \n')
                for id, xco, yco, zco, L in zip(ID, x, y, z, Labs):
                    f.write(
                        str(id) + '   ' + str(xco) + '   ' + str(yco) + '   ' + str(zco) + '   ' + str(L) + '\n')

        def writeLabs(file, ID, Labs):

            with open(file, 'w') as f:
                f.write('# ID    Labs(W) \n')
                for id, L in zip(ID, Labs):
                    f.write(str(id) + '   ' + str(L) + '\n')

    # -----------------------------------------------------------------

    def load_absorption(self):

        """
        This function ...
        :return:
        """

        # _ds_abs.dat

        # Bolometric absorbed luminosituy for all dust cells with nonzero absorption
        # column 1: dust cell index
        # column 2: x coordinate of cell center (pc)
        # column 3: y coordinate of cell center (pc)
        # column 4: z coordinate of cell center (pc)
        # column 5: L_abs (W)

        # Loop over the different contributions
        for contribution in self.contributions:

            # Determine the path to the output directory of the simulation
            output_path = self.output_paths[contribution]

            # Determine the path to the ISFR data file
            isrf_path = fs.join(output_path, self.galaxy_name + "_ds_abs.dat")

            # Load the ISRF file
            cell, x, y, labs = np.loadtxt(isrf_path, usecols=(0,1,2,3,4), unpack=True)



    # -----------------------------------------------------------------

    def integrate_over_wavelengths(self, jlambdas):

        """
        This function ...
        :param jlambdas:
        :return:
        """

        lum_cells = []

        # L = sum_lambda (j_lambda * d_lambda)

        # Get the wavelength deltas
        deltas = self.wavelength_grid["Delta"]
        deltas = tables.column_as_array(deltas, unit="m") # deltas in meter

        # Loop over all dust cells
        for i in range(jlambdas[0].size):

            # Put the Jlambda values for this cell in an array
            jlambdas_cell = np.array([jlambdas[k][i] for k in range(len(jlambdas))])

            # Calculate the luminosity
            #lum = np.sum(jlambdas_cell * MjySr_to_LsunMicron * deltas)
            lum = np.sum(jlambdas_cell * deltas) # if deltas are in meter, this is in W/m3/sr * m -> W/m2/sr
            #lum = np.sum(jlambdas_cell * deltas * ...)

            # Add the luminosity to the list
            lum_cells.append(lum)

        # Return the list of luminosities
        return lum_cells

    # -----------------------------------------------------------------

    def merge(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the PTS modeling/analysis directory
        pts_modeling_analysis_path = fs.join(inspection.pts_package_dir, "modeling", "analysis")

        # Determine the path to the R script that does the merge
        #script_path = fs.join(pts_modeling_analysis_path, "merge_sets.R")
        #subprocess.call(["Rscript", script_path])

        # Define the lines that have to be executed by R to merge the absorption data
        lines = []
        lines.append("library('dplyr')")
        lines.append("path = './")
        lines.append("props < - read.table(paste0(path, 'M31_212isrf_new_ds_cellprops.dat'))")
        lines.append("IDs < - c(0:(length(props[[1]]) - 1))")
        lines.append("props < - mutate(props, ID=IDs)")

        lines.append("tot < - read.table(paste0(path, 'Labs_tot.dat'))")
        lines.append("old < - read.table(paste0(path, 'Labs_old.dat'))")
        lines.append("yng < - read.table(paste0(path, 'Labs_yng.dat'))")
        lines.append("new < - read.table(paste0(path, 'Labs_new.dat'))")

        lines.append("propstot < - merge(props, tot, by.x = 'ID', by.y = 'V1')")
        lines.append("print('Merged total properties...')")

        lines.append("oldyng < - merge(old, yng, by.x = 'V1', by.y = 'V1')")
        lines.append("oldyngnew < - merge(oldyng, new, by.x = 'V1', by.y = 'V1')")
        lines.append("print('Merged component properties...')")

        lines.append("final < - merge(propstot, oldyngnew, by.x = 'ID', by.y = 'V1')")
        lines.append("names(final) < - c('ID', 'volume', 'density', 'massFraction', 'odepth', 'density_new', 'x', 'y', 'z', 'Ltot', 'Lold', 'Lyng', 'Lnew')")
        lines.append("print('Final set created...')")

        lines.append("write.table(final, paste0(path, 'total212Labs.dat'), row.names = F)")

        # Execute all lines consecutively in R
        for line in lines: ro.r(line)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def seba(self):

        outpath = "modelChecks/iteration5_J14/"
        inpath = "SKIRTOutput/iteration5_J14/"
        inSED = "M31_212full_i77.5_sed.dat"

        Lsun = 3.846e26  # Watts

        # load SEDs
        input = np.loadtxt(inpath + inSED)
        wavelengths = input[:, 0]

        # Load the widths of the wavelength bins. Crucial for integration!
        delta_wls = np.loadtxt("SKIRTOutput/iteration5_J14/M31_reference_wavelengths512.dat", usecols=(1,))

        # only consider wavelengths longwards of 10 micron. For speed and memory
        startwl = next(wl[0] for wl in enumerate(wavelengths) if wl[1] > 10.)
        coldstartwl = next(wl[0] for wl in enumerate(wavelengths) if wl[1] > 100.)

        # produce wavelength ranges for all, warm and cold dust.
        dustwls = wavelengths[startwl:]
        warmwls = wavelengths[startwl:coldstartwl - 1]
        coldwls = wavelengths[coldstartwl:]
        delta_wls = delta_wls[startwl:]

        # Compute global heating fracions

        flux_all = input[startwl:, 1]

        input = np.loadtxt(inpath + inSED.replace('_i', '_old_i'))
        flux_old = input[startwl:, 1]

        input = np.loadtxt(inpath + inSED.replace('_i', '_young_i'))
        flux_young = input[startwl:, 1]

        Fold = 100. * (0.5 * flux_old + 0.5 * (flux_all - flux_young)) / flux_all
        Fyoung = 100. * (0.5 * flux_young + 0.5 * (flux_all - flux_old)) / flux_all

        Fold_alternative1 = 100. * flux_old / (flux_old + flux_young)
        Fyoung_alternative1 = 100. * flux_young / (flux_old + flux_young)

        Fold_alternative2 = 100. * flux_old / flux_all
        Fyoung_alternative2 = 100. * flux_young / flux_all

        Fold_alternative3 = 100. * np.sqrt(flux_old * (flux_all - flux_young)) / flux_all
        Fyoung_alternative3 = 100. * np.sqrt(flux_young * (flux_all - flux_old)) / flux_all

        Fold_alternative4 = 100. * (0.5 * flux_old + 0.5 * (flux_all - flux_young)) / (flux_old + flux_young)
        Fyoung_alternative4 = 100. * (0.5 * flux_young + 0.5 * (flux_all - flux_old)) / (flux_old + flux_young)

        JyToLsun = 1.e-26 * 4 * np.pi * (0.785e6 * 3.086e+16) ** 2 * 3.e14 / (dustwls ** 2) / Lsun  # in Lsun/micron

        totFlux = 0
        totFlux_young = 0
        totFlux_old = 0
        for i in range(len(flux_all) - 1):
            totFlux += delta_wls[i] * flux_all[i] * JyToLsun[i]
            totFlux_young += delta_wls[i] * flux_young[i] * JyToLsun[i]
            totFlux_old += delta_wls[i] * flux_old[i] * JyToLsun[i]

        print('Total heating from old stars: ', totFlux_old / totFlux)
        print('Total heating from young stars: ', totFlux_young / totFlux)

        plt.figure(figsize=(7, 5))
        plt.ylabel('$F^\prime_{\lambda,\mathrm{unev.}}$ [$\%$]', fontsize=20)
        plt.xlabel('$\lambda/\mu\mathrm{m}$', fontsize=20)
        plt.xlim(10., 1.e3)
        plt.ylim(0., 60.)
        plt.xscale('log')
        plt.tick_params(labelsize=20)
        # plt.subplots_adjust(bottom=0.1)
        # plt.plot(dustwls,Fold, 'r-', label="old stars")
        plt.plot(dustwls, Fyoung, 'k-', label="Young SPs")
        # plt.plot(dustwls,Fyoung_alternative1, 'r-', label="alt 1")
        # plt.plot(dustwls,Fyoung_alternative2, 'g-', label="alt 2")
        # plt.plot(dustwls,Fyoung_alternative3, 'c-', label="alt 3")
        plt.plot(dustwls, Fyoung_alternative4, 'k-', label="alt 4")
        plt.fill_between(dustwls, Fyoung, Fyoung_alternative4, color='grey', alpha='0.5')

        plt.tight_layout()
        # plt.legend(loc='upper left',numpoints=1,markerscale=1.5,fontsize=14)
        plt.savefig(outpath + inSED.replace('sed.dat', 'heating.pdf'), format='pdf')

        # plt.show()
        plt.close()

        # Make heating map
        inCube = inSED.replace('sed.dat', 'total.fits')
        makeHeatMap(inpath, outpath, inCube, startwl, dustwls, delta_wls)
        # makeWarmHeatMap(inpath,outpath,inCube,startwl,coldstartwl-1, warmwls)
        # makeColdHeatMap(inpath,outpath,inCube,coldstartwl, coldwls)

# -----------------------------------------------------------------

def makeHeatMap(inpath,outpath,inCube,startwl,dustwls, delta_wls):

    cube = pyfits.open(inpath+inCube)
    cube_all = cube[0].data[startwl:,0:,0:]
    hdr_all  = cube[0].header

    cube = pyfits.open(inpath+inCube.replace('_i','_old_i'))
    cube_old = cube[0].data[startwl:,0:,0:]

    cube = pyfits.open(inpath+inCube.replace('_i','_young_i'))
    cube_young = cube[0].data[startwl:,0:,0:]


    Fold   = 100. * (0.5*cube_old + 0.5*(cube_all-cube_young)) / cube_all
    hdu = pyfits.PrimaryHDU( Fold,hdr_all)
    hdu.writeto(outpath+"heatingFold.fits",clobber=True)

    Fyoung = 100. * (0.5*cube_young + 0.5*(cube_all-cube_old)) / cube_all
    hdu = pyfits.PrimaryHDU( Fyoung,hdr_all)
    hdu.writeto(outpath+"heatingFyoung.fits",clobber=True)


    pixelTot = integratePixelSEDs(cube_all,     dustwls, delta_wls)
    pixelOld = integratePixelSEDs(cube_old,     dustwls, delta_wls)
    pixelYoung = integratePixelSEDs(cube_young, dustwls, delta_wls)

    # Get header with appropriate WCS info
    im36 = pyfits.open("SKIRTinput/new3.6MJySr.fits")
    hdr_wcs = im36[0].header

    hdu = pyfits.PrimaryHDU(pixelTot,hdr_wcs)
    hdu.writeto(outpath+"Ldust_tot.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelOld,hdr_wcs)
    hdu.writeto(outpath+"Ldust_old.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelYoung,hdr_wcs)
    hdu.writeto(outpath+"Ldust_young.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelOld/pixelTot,hdr_wcs)
    hdu.writeto(outpath+"heatingTotOld.fits",clobber=True)

    hdu = pyfits.PrimaryHDU(pixelYoung/pixelTot,hdr_wcs)
    hdu.writeto(outpath+"heatingTotYoung.fits",clobber=True)


# OLD AND INCORRECT?
#    tot_all = 100.* (dustwls[len(dustwls)-1] - dustwls[0])
#
#    tot_old   = integrateHeating(Fold,dustwls) / tot_all
#    hdu = pyfits.PrimaryHDU(tot_old,hdr_wcs)
#    hdu.writeto(outpath+"heatingTotOld.fits",clobber=True)
#
#    tot_young = integrateHeating(Fyoung,dustwls) / tot_all
#    hdu = pyfits.PrimaryHDU(tot_young,hdr_wcs)
#    hdu.writeto(outpath+"heatingTotYoung.fits",clobber=True)

def makeWarmHeatMap(inpath,outpath,inCube,startwl,stopwl,warmwls):

    cube = pyfits.open(inpath+inCube)
    cube_all = cube[0].data[startwl:stopwl,0:,0:]
    hdr_all  = cube[0].header

    cube = pyfits.open(inpath+inCube.replace('_i','_old_i'))
    cube_old = cube[0].data[startwl:stopwl,0:,0:]

    cube = pyfits.open(inpath+inCube.replace('_i','_young_i'))
    cube_young = cube[0].data[startwl:stopwl,0:,0:]

    Fold   = 100. * (0.5*cube_old + 0.5*(cube_all-cube_young)) / cube_all
    Fyoung = 100. * (0.5*cube_young + 0.5*(cube_all-cube_old)) / cube_all

    tot_all = 100.* (warmwls[len(warmwls)-1] - warmwls[0])

    # Get header with appropriate WCS info
    im36 = pyfits.open("SKIRTinput/new3.6MJySr.fits")
    hdr_wcs = im36[0].header

    tot_old   = integrateHeating(Fold,warmwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_old,hdr_wcs)
    hdu.writeto(outpath+"heatingTotWarmOld.fits",clobber=True)

    tot_young = integrateHeating(Fyoung,warmwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_young,hdr_wcs)
    hdu.writeto(outpath+"heatingTotWarmYoung.fits",clobber=True)

def makeColdHeatMap(inpath,outpath,inCube,startwl,coldwls):

    cube = pyfits.open(inpath+inCube)
    cube_all = cube[0].data[startwl:,0:,0:]
    hdr_all  = cube[0].header

    cube = pyfits.open(inpath+inCube.replace('_i','_old_i'))
    cube_old = cube[0].data[startwl:,0:,0:]

    cube = pyfits.open(inpath+inCube.replace('_i','_young_i'))
    cube_young = cube[0].data[startwl:,0:,0:]


    Fold   = 100. * (0.5*cube_old + 0.5*(cube_all-cube_young)) / cube_all
    Fyoung = 100. * (0.5*cube_young + 0.5*(cube_all-cube_old)) / cube_all

    tot_all = 100.* (coldwls[len(coldwls)-1] - coldwls[0])

    # Get header with appropriate WCS info
    im36 = pyfits.open("SKIRTinput/new3.6MJySr.fits")
    hdr_wcs = im36[0].header

    tot_old   = integrateHeating(Fold,coldwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_old,hdr_wcs)
    hdu.writeto(outpath+"heatingTotColdOld.fits",clobber=True)

    tot_young = integrateHeating(Fyoung,coldwls) / tot_all
    hdu = pyfits.PrimaryHDU(tot_young,hdr_wcs)
    hdu.writeto(outpath+"heatingTotColdYoung.fits",clobber=True)


def integratePixelSEDs(cube, wls, dwls):


    Lsun = 3.846e26 # Watts
    MjySr_to_LsunMicron = 1.e6 * (36./206264.806247)**2 * 1.e-26 * 4*np.pi*(0.785e6*3.086e+16)**2 * 3.e14/(wls**2) / Lsun

    xaxis = len(cube[0,0,0:])
    yaxis = len(cube[0,0:,0])
    zaxis = len(cube[0:,0,0])

    slice = np.zeros((yaxis,xaxis))
    for i in range(0,yaxis):
        for j in range(0,xaxis):
            sed = cube[0:,i,j]
            slice[i,j] = np.sum(sed * MjySr_to_LsunMicron * dwls)

    return slice


def integrateHeating(cube,dustwls):

    xaxis = len(cube[0,0,0:])
    yaxis = len(cube[0,0:,0])
    zaxis = len(cube[0:,0,0])

    slice = np.zeros((yaxis,xaxis))
    for i in range(0,yaxis):
        for j in range(0,xaxis):
            sed = cube[0:,i,j]
            slice[i,j] = integrate.simps(sed,dustwls)

    return slice

# -----------------------------------------------------------------
