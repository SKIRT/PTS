#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_correlations Plot correlations between different properties of a galaxy table.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.table import SmartTable
from pts.core.tools import filesystem as fs
from pts.core.tools import sequences
from pts.magic.tools import plotting

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Required
definition.add_required("filename", "file_path", "data table filename")
definition.add_optional("output", "directory_path", "output directory path")

# Plotting format
definition.add_optional("format", "string", "plotting format", "pdf")

# Replot
definition.add_optional("replot", "string_list", "remake plots with these variables")

# Get configuration
config = parse_arguments("plot_correlations", definition)

# -----------------------------------------------------------------

# Determine full path
filepath = fs.absolute_or_in_cwd(config.filename)

# Load the table
table = SmartTable.from_file(filepath)
galaxy_names = list(table["name"])

# -----------------------------------------------------------------

lums = ["galex_fuv", "galex_nuv", "sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z", "2mass_j", "2mass_h", "2mass_k", "wise_w1", "wise_w2", "wise_w3", "wise_w4", "iras12", "iras25", "iras60", "iras100", "i1", "i2", "i3", "i4", "mips24", "mips70", "mips160", "pblue", "pgreen", "pred", "psw", "pmw", "plw", "hfi_350", "hfi_550", "hfi_850", "hfi_1380", "hfi_2100", "hfi_3000"]
colours = ["fuv_nuv", "fuv_h", "fuv_j", "fuv_k", "fuv_u", "fuv_g", "fuv_r", "fuv_i", "fuv_z", "fuv_mu3", "fuv_mu4", "nuv_h", "nuv_j", "nuv_k", "nuv_u", "nuv_g", "nuv_r", "nuv_i", "nuv_z", "nuv_mu3", "nuv_mu4", "mu25_mu70", "mu25_mu60", "mu60_mu100", "mu70_mu100", "mu100_mu160", "mu160_mu250", "mu250_mu350", "mu350_mu500", "mu350_mu550", "mu500_mu850", "mu550_mu850", "mu850_mu1380", "mu1380_mu2100", "mu2100_mu3000"]

# -----------------------------------------------------------------

filternames_for_colours = dict()
filternames_for_colours["fuv"] = ["galex_fuv"]
filternames_for_colours["nuv"] = ["galex_nuv"]
filternames_for_colours["h"] = ["2mass_h"]
filternames_for_colours["j"] = ["2mass_j"]
filternames_for_colours["k"] = ["2mass_k"]
filternames_for_colours["u"] = ["sdss_u"]
filternames_for_colours["g"] = ["sdss_g"]
filternames_for_colours["r"] = ["sdss_r"]
filternames_for_colours["i"] = ["sdss_i"]
filternames_for_colours["z"] = ["sdss_z"]
filternames_for_colours["mu3"] = ["wise_w1", "i1"]
filternames_for_colours["mu4"] = ["wise_w2", "i2"]
filternames_for_colours["mu25"] = ["iras25", "mips24"]
filternames_for_colours["mu70"] = ["pblue", "mips70"]
filternames_for_colours["mu60"] = ["iras60"]
filternames_for_colours["mu100"] = ["iras100", "pgreen"]
filternames_for_colours["mu160"] = ["pred", "mips160"]
filternames_for_colours["mu250"] = ["psw"]
filternames_for_colours["mu350"] = ["pmw", "hfi_350"]
filternames_for_colours["mu500"] = ["plw"]
filternames_for_colours["mu550"] = ["hfi_550"]
filternames_for_colours["mu850"] = ["hfi_850"]
filternames_for_colours["mu1380"] = ["hfi_1380"]
filternames_for_colours["mu2100"] = ["hfi_2100"]
filternames_for_colours["mu3000"] = ["hfi_3000"]

# -----------------------------------------------------------------

ignore_columns = ["Name", "Common name", "RA", "DEC", "Inclination", "Ellipticity", "Redshift", "Names", "Distance", "Position angle", "Velocity", "D25", "Effective radius", "Type"]
column_names = [column_name for column_name in table.colnames if column_name not in ignore_columns and not column_name.startswith("has_") and not column_name.endswith("error")]

# -----------------------------------------------------------------

linear_scales = ["stage"]

#print(column_names)
#print(len(column_names))

# -----------------------------------------------------------------

ncombinations = 0
for column_a, column_b in sequences.combinations(column_names, 2):

    a_is_colour = column_a in colours
    b_is_colour = column_b in colours

    if a_is_colour and b_is_colour: continue

    a_is_lum = column_a in lums
    b_is_lum = column_b in lums

    if a_is_lum and b_is_lum: continue

    if a_is_lum and b_is_colour:

        b_part1, b_part2 = column_b.split("_")
        #b_filters = get_filters_for_colour(column_b)
        b_filter_names = filternames_for_colours[b_part1] + filternames_for_colours[b_part2]
        if column_a in b_filter_names: continue

    elif a_is_colour and b_is_lum:

        a_part1, a_part2 = column_a.split("_")
        #a_filters = get_filters_for_colour(column_a)
        a_filter_names = filternames_for_colours[a_part1] + filternames_for_colours[a_part2]
        if column_b in a_filter_names: continue

    #print(column_a, column_b)

    # Set filename
    filename = column_a + "__" + column_b
    output_path = fs.join(config.output_path(), filename + "." + config.format)

    # Check
    if fs.is_file(output_path):
        if config.replot is not None and (column_a in config.replot or column_b in config.replot): fs.remove_file(output_path)
        else: continue

    # Logscales
    if a_is_colour or column_a in linear_scales: xlog = False
    else: xlog = True
    if b_is_colour or column_b in linear_scales: ylog = False
    else: ylog = True

    # Plot
    plotting.plot_table(filepath, column_a, column_b, output_path, x_log=xlog, y_log=ylog)

    ncombinations += 1

print(ncombinations)

# topcat -stilts plot2plane \
#    xpix=1156 ypix=593 \
#    xlog=true ylog=true xlabel=mips70 ylabel=sfr xcrowd=0.9998301109057076 \
#     ycrowd=0.9998301109057076 \
#    xmin=7.0E26 xmax=1.0498369777416409E35 ymin=2.0E-9 ymax=30.260115332352253 \
#    auxmin=1.0278138959234115 auxmax=4.879091762160679 \
#    auxvisible=true auxlabel=nuv_g auxcrowd=0.9998301109057076 \
#    legend=false \
#    in=/Users/samverstocken/MODELING/galaxies2.dat ifmt=ASCII x=mips70 y=sfr \
#    layer_1=Mark \
#       aux_1=nuv_g \
#       shading_1=aux size_1=5 \
#    layer_2=LinearFit out=plot.pdf

# -----------------------------------------------------------------
