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
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import rpy2.robjects as ro

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ....core.tools import tables, introspection
from ....core.simulation.table import SkirtTable
from ....core.basics.distribution import Distribution, Distribution2D
from ....core.simulation.wavelengthgrid import WavelengthGrid

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

        # -- Attributes --

        # The wavelength grid used for the simulations
        self.wavelength_grid = None

        # The number of wavelengths
        self.number_of_wavelengths = None

        # The table with the cell properties
        self.cell_properties = None

        # The table with the absorbed luminosities
        self.absorptions = None

        # The mask of cells for which the total absorbed luminosity is zero
        self.zero_absorption = None

        # The heating fraction of the unevolved stellar population for each dust cell
        self.heating_fractions = None

        # The distribution of heating fractions
        self.distribution = None

        # The 2D distribution of heating fractions
        self.radial_distribution = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the wavelength grid
        self.load_wavelength_grid()

        # 3. Load the cell properties
        self.load_cell_properties()

        # 4. Load the absorption data
        self.load_absorption()

        # 5. Calculate the heating fraction of the unevolved stellar population
        self.calculate_heating_unevolved()

        # 6. Calculate the distribution of the heating fraction of the unevolved stellar population
        self.calculate_distribution()

        # 7. Calculate the distribution of the heating fraction of the unevolved stellar population as a function of radius
        self.calculate_radial_distribution()

        # 8. Writing
        self.write()

        # 9. Plotting
        self.plot()

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
        self.wavelength_grid = WavelengthGrid.from_skirt_output(wavelengths_path)

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

        # Load the properties table
        self.cell_properties = SkirtTable.from_file(properties_path)

    # -----------------------------------------------------------------

    def load_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the absorption data ...")

        # Header:
        # Bolometric absorbed luminosity for all dust cells
        # column 1: x coordinate of cell center (pc)
        # column 2: y coordinate of cell center (pc)
        # column 3: z coordinate of cell center (pc)
        # column 4: Absorbed bolometric luminosity (W)

        contribution_tables = dict()

        # Loop over the different contributions
        for contribution in self.contributions:

            # Skip the simulation of the total unevolved (young + ionizing) stellar population
            if contribution == "unevolved": continue

            # Debugging
            log.debug("Loading the SKIRT absorption table for the simulation of the " + contribution + " stellar population ...")

            # Determine the path to the output directory of the simulation
            output_path = self.output_paths[contribution]

            # Determine the path to the absorption data file
            absorption_path = fs.join(output_path, self.galaxy_name + "_ds_abs.dat")

            # Load the absorption table for this contribution
            table = SkirtTable.from_file(absorption_path)

            # Add the table
            contribution_tables[contribution] = table

        do_checks = False
        if do_checks:

            # Debugging
            log.debug("Checking whether the tables are consistent ...")

            # Check whether the table lengths match
            table_lengths = [len(table) for table in contribution_tables.values()]
            if not all(table_lengths[0] == length for length in table_lengths): raise ValueError("Absorption tables have different sizes")

            # Check whether the X coordinates of the cells match
            if not tables.equal_columns([table["X coordinate of cell center"] for table in contribution_tables.values()]):
                raise ValueError("Columns of X coordinates of cell centers do not match between the different contributions")

            # Check whether the Y coordinates of the cells match
            if not tables.equal_columns([table["Y coordinate of cell center"] for table in contribution_tables.values()]):
                raise ValueError("Columns of Y coordinates of cell centers do not match between the different contributions")

            # Check whether the Z coordinates of the cells match
            if not tables.equal_columns([table["Z coordinate of cell center"] for table in contribution_tables.values()]):
                raise ValueError("Columns of Z coordinates of cell centers do not match between the different contributions")

        # Debugging
        log.debug("Creating the absorption table ...")

        # Create the columns for the absorption table
        #names = ["X coordinate of cell center", "Y coordinate of cell center", "Z coordinate of cell center"]
        #data = []
        #data.append(contribution_tables[contribution_tables.keys()[0]]["X coordinate of cell center"])
        #data.append(contribution_tables[contribution_tables.keys()[0]]["Y coordinate of cell center"])
        #data.append(contribution_tables[contribution_tables.keys()[0]]["Z coordinate of cell center"])

        # Loop over the tables of the different contributions, add the absorption luminosity columns
        #for contribution in contribution_tables:

            #names.append("Absorbed bolometric luminosity of the " + contribution + " stellar population")
            #data.append(contribution_tables[contribution]["Absorbed bolometric luminosity"])

        # Create the absorption table
        #self.absorptions = tables.new(data, names, copy=False)

        self.absorptions = tables.new()

        x_coords = contribution_tables[contribution_tables.keys()[0]]["X coordinate of cell center"]
        y_coords = contribution_tables[contribution_tables.keys()[0]]["Y coordinate of cell center"]
        z_coords = contribution_tables[contribution_tables.keys()[0]]["Z coordinate of cell center"]
        self.absorptions.add_columns([x_coords, y_coords, z_coords], copy=False)

        for contribution in contribution_tables:

            self.absorptions.add_columns([contribution_tables[contribution]["Absorbed bolometric luminosity"]], copy=False)
            self.absorptions.rename_column("Absorbed bolometric luminosity", "Absorbed bolometric luminosity of the " + contribution + " stellar population")

        # Create a mask of cells with zero absorption
        self.zero_absorption = self.absorptions["Absorbed bolometric luminosity of the total stellar population"] == 0.

    # -----------------------------------------------------------------

    def calculate_heating_unevolved(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the heating fraction of the unevolved stellar population ...")

        # SEBA:

        #totISRFfile = "total212Labs.dat"
        #Lsun = 3.846e26  # Watts

        ## Total energy absorbed in the new component
        ## Derived from ../../SKIRTrun/models/testHeating/MappingsHeating/plotSEDs.py
        #Lnew = 176495776.676  # in Lsun

        #print ''reading data... '
        #input = np.loadtxt(path + totISRFfile, skiprows=1)
        #ID = input[:, 0]
        #volume = input[:, 1]
        #density = input[:, 2]
        #massFrac = input[:, 3]
        #density_new = input[:, 5]
        #x = input[:, 6]
        #y = input[:, 7]
        #tot = input[:, 9] / Lsun
        #old = input[:, 10] / Lsun
        #yng = input[:, 11] / Lsun
        #new = input[:, 12] / Lsun

        #energy_new = volume * density_new * Lnew
        #F_abs_yng = (yng + new + energy_new) / (old + yng + new + energy_new)

        absorptions_unevolved_diffuse = self.absorptions["Absorbed bolometric luminosity of the young stellar population"] + self.absorptions["Absorbed bolometric luminosity of the ionizing stellar population"]
        #absorptions_ionizing_internal = None # TODO !!

        absorptions_total = self.absorptions["Absorbed bolometric luminosity of the total stellar population"]
        #absorptions_total = absorptions_unevolved_diffuse + absorptions_ionizing_internal + absorptions_evolved # TODO !!

        # Calculate the heating fraction of the unevolved stellar population in each dust cell
        self.heating_fractions = absorptions_unevolved_diffuse / absorptions_total



        self.mask = self.heating_fractions.mask  # is basically the same mask as self.zero_absorption because during the division (see calculate_heating_unevolved) Astropy sets invalid values as masked
        # nan_inf_mask = np.isnan(self.heating_fractions) + np.isinf(self.heating_fractions) # no nans or infs because Astropy sets them as MaskedConstants during the division
        greater_than_one_mask = self.heating_fractions > 1.0

        # Debugging
        log.debug(str(np.sum(greater_than_one_mask)) + " pixels have a heating fraction greater than unity")

        self.mask += greater_than_one_mask
        self.heating_fractions_compressed = np.ma.MaskedArray(self.heating_fractions, self.mask).compressed()
        self.weights_compressed = np.ma.MaskedArray(self.cell_properties["Mass fraction"], self.mask).compressed()

    # -----------------------------------------------------------------

    def calculate_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population ...")

        # Generate the distribution
        self.distribution = Distribution.from_values(self.heating_fractions_compressed, bins=20, weights=self.weights_compressed)

    # -----------------------------------------------------------------

    def calculate_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the radial distribution of heating fractions of the unevolved stellar population ...")

        # Calculate the radius for each dust cell
        x_coords = self.absorptions["X coordinate of cell center"]
        y_coords = self.absorptions["Y coordinate of cell center"]
        radii = np.sqrt(x_coords**2 + y_coords**2)

        radii_compressed = np.ma.MaskedArray(radii, self.mask).compressed()

        # Generate the radial distribution
        self.radial_distribution = Distribution2D.from_values(radii_compressed, self.heating_fractions_compressed, weights=self.weights_compressed, x_name="radius (pc)", y_name="Heating fraction of unevolved stars")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the distribution of heating fractions
        self.write_distribution()

        # Write the radial distribution of heating fractions
        self.write_radial_distribution()

    # -----------------------------------------------------------------

    def write_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distribution of heating fractions of the unevolved stellar population ...")

        # Determine the path to the distribution file
        path = fs.join(self.analysis_heating_path, "distribution.dat")

        # Save the distribution
        self.distribution.save(path)

    # -----------------------------------------------------------------

    def write_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the radial distribution of heating fractions of the unevolved stellar population ...")

        # Determine the path to the radial distribution file
        path = fs.join(self.analysis_heating_path, "radial_distribution.dat")

        # Save the radial distribution
        self.radial_distribution.save(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the distribution of heating fractions
        self.plot_distribution()

        # Plot the radial distribution of heating fractions
        self.plot_radial_distribution()

        # Plot a map of the heating fraction for a face-on view of the galaxy
        self.plot_map()

    # -----------------------------------------------------------------

    def plot_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a histogram of the distribution of heating fractions of the unevolved stellar population ...")

        # Determine the path to the plot file
        path = fs.join(self.analysis_heating_path, "distribution.pdf")

        # Create the plot file
        self.distribution.plot(title="Distribution of the heating fraction of the unevolved stellar population", path=path)

    # -----------------------------------------------------------------

    def plot_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a 2D histogram of the radial distribution of the heating fractions of the unevolved stellar population ...")

        # Determine the path to the plot file
        path = fs.join(self.analysis_heating_path, "radial_distribution.pdf")

        # Create the plot file
        self.radial_distribution.plot(title="Radial distribution of the heating fraction of the unevolved stellar population", path=path)

    # -----------------------------------------------------------------

    def plot_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a map of the heating fraction of the unevolved stellar population for a face-on view of the galaxy ...")

        # Determine the path to the plot file
        path = fs.join(self.analysis_heating_path, "map.pdf")

        plt.figure()

        x = np.ma.MaskedArray(self.absorptions["X coordinate of cell center"], mask=self.mask).compressed()
        y = np.ma.MaskedArray(self.absorptions["Y coordinate of cell center"], mask=self.mask).compressed()
        z = self.heating_fractions_compressed

        #plt.pcolormesh(x, y, z, cmap='RdBu', vmin=0.0, vmax=1.0)

        from matplotlib import mlab

        x_ticks = x
        y_ticks = y

        z_grid = mlab.griddata(x, y, z, x, y)

        from matplotlib.backends import backend_agg as agg
        from matplotlib import cm

        # plot
        # fig = Figure()  # create the figure
        fig = plt.figure()
        agg.FigureCanvasAgg(fig)  # attach the rasterizer
        ax = fig.add_subplot(1, 1, 1)  # make axes to plot on
        ax.set_title("Interpolated Contour Plot of Experimental Data")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        cmap = cm.get_cmap("hot")  # get the "hot" color map
        contourset = ax.contourf(x_ticks, y_ticks, z_grid, 10, cmap=cmap)

        plt.savefig(path)
        plt.close()

    # -----------------------------------------------------------------

    def load_absorption_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the absorption data ...")

        # Bolometric absorbed luminosituy for all dust cells with nonzero absorption
        # column 1: dust cell index
        # column 2: x coordinate of cell center (pc)
        # column 3: y coordinate of cell center (pc)
        # column 4: z coordinate of cell center (pc)
        # column 5: L_abs (W)

        data = dict()

        # Loop over the different contributions
        for contribution in self.contributions:

            # Debugging
            log.debug("Loading the SKIRT absorption table for the simulation of the " + contribution + " stellar population ...")

            # Determine the path to the output directory of the simulation
            output_path = self.output_paths[contribution]

            # Determine the path to the absorption data file
            absorption_path = fs.join(output_path, self.galaxy_name + "_ds_abs.dat")

            # Load the absorption table for this contribution
            table = SkirtTable.from_file(absorption_path)

            # Add the table
            data[contribution] = table

        # Create the absorption table which contains the information for the different contributions

        # Debugging
        log.info("Merging the individual absorption tables into one (this can take a while) ...")

        current_rows = dict()
        for contribution in self.contributions: current_rows[contribution] = 0

        x_column = []
        y_column = []
        z_column = []
        abs_columns = defaultdict(list)

        # Loop over the cells
        number_of_cells = len(self.cell_properties)
        for index in range(number_of_cells):

            # Debugging
            log.debug("Cell " + str(index+1) + " of " + str(number_of_cells) + " (" + str((index+1)/number_of_cells*100.) + "%) ...")

            x = None
            y = None
            z = None
            abs = dict()

            # Loop over the absorption tables of the different contributions, find the cell index
            for contribution in self.contributions:

                table = data[contribution]
                current_row = current_rows[contribution]

                if table["Dust cell index"][current_row] == index:

                    x = table["X coordinate of cell center"][current_row]
                    y = table["Y coordinate of cell center"][current_row]
                    z = table["Z coordinate of cell center"][current_row]
                    labs = table["L_abs"][current_row]
                    abs[contribution] = labs

                    # Use the next row to test against the next cell index
                    current_rows[contribution] += 1

            x_column.append(x)
            y_column.append(y)
            z_column.append(z)

            for contribution in self.contributions:
                abs_columns[contribution].append(abs[contribution] if contribution in abs else None)

        # Initialize the column data and names
        data = [x_column, y_column, z_column]
        names = ["X coordinate", "Y coordinate", "Z coordinate"]
        for contribution in self.contributions:
            data.append(abs_columns[contribution])
            names.append("Absorption for " + contribution + " stellar population")

        # Create the final absorption table
        self.absorptions = tables.new(data, names)

    # -----------------------------------------------------------------

    def load_absorption_in_r(self):

        """
        This function ...
        :return:
        """

        properties_path = fs.join(self.output_paths["total"], self.galaxy_name + "_ds_cellprops.dat")
        tot_path = fs.join(self.output_paths["total"], self.galaxy_name + "_ds_abs.dat")
        old_path = fs.join(self.output_paths["old"], self.galaxy_name + "_ds_abs.dat")
        young_path = fs.join(self.output_paths["young"], self.galaxy_name + "_ds_abs.dat")
        ionizing_path = fs.join(self.output_paths["ionizing"], self.galaxy_name + "_ds_abs.dat")
        self.merge_in_r(properties_path, tot_path, old_path, young_path, ionizing_path)

    # -----------------------------------------------------------------

    def merge_in_r(self, properties_path, tot_path, old_path, young_path, ionizing_path):

        """
        This function ...
        :return:
        """

        temp_path = fs.join(self.analysis_heating_path, "abs_temp.dat")

        # Define the lines that have to be executed by R to merge the absorption data
        lines = []
        lines.append("library('dplyr')")
        lines.append("props <- read.table('" + properties_path + "')")
        lines.append("IDs <- c(0:(length(props[[1]]) - 1))")
        lines.append("props <- mutate(props, ID=IDs)")

        lines.append("tot <- read.table('" + tot_path + "')")
        lines.append("old <- read.table('" + old_path + "')")
        lines.append("yng <- read.table('" + young_path + "')")
        lines.append("new <- read.table('" + ionizing_path + "')")

        lines.append("propstot <- merge(props, tot, by.x = 'ID', by.y = 'V1')")
        lines.append("print('Merged total properties...')")

        lines.append("oldyng <- merge(old, yng, by.x = 'V1', by.y = 'V1')")
        lines.append("oldyngnew <- merge(oldyng, new, by.x = 'V1', by.y = 'V1')")
        lines.append("print('Merged component properties...')")

        lines.append("final <- merge(propstot, oldyngnew, by.x = 'ID', by.y = 'V1')")
        lines.append("names(final) <- c('ID', 'volume', 'density', 'massFraction', 'odepth', 'density_new', 'x', 'y', 'z', 'Ltot', 'Lold', 'Lyng', 'Lnew')")
        lines.append("print('Final set created...')")

        lines.append("write.table(final, '" + temp_path + "', row.names = F)")

        # Execute all lines consecutively in R
        for line in lines:
            #print(line)
            ro.r(line)

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
        deltas = self.wavelength_grid.deltas(asarray=True, unit="m") # deltas in meter

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
