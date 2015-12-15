#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.fitskirt_results
# Create a file containing results for a polychromatic FitSKIRT fit.
# The only input is the fski files used to create the final results.
# This script creates both the separate result for each fit as well 
# as the statistical best fit and error for all best_simulation 
# files in current folder or subfolders.
#
# REQUIREMENTS: 
# - .fski used to fit the models. Only results of this .fski (or copies) will be included.
#
# EXAMPLE:>pts fitskirt_results.py ranges.fski
##

# ------------------------------------------------------------------

# Import standard modules
import sys
import glob
import os
import math

# ------------------------------------------------------------------

# fski files used to read the labels and units
fski_path=str(sys.argv[1])
fskiname = fski_path.split('.fski')[0].split('/')[len(fski_path.split('.fski')[0].split('/'))-1]

# fetching labels, units of the free parameters and the skifile
skifile =''
labels, Units, wavelengths, matches, solutions = [], [], [], [], []
fski = open(fski_path,'r')

for line in fski:
	if '<AdjustableSkirtSimulation' in line: skifile=line.split('"')[len(line.split('"'))-2]
	if '<ParameterRange label=' in line:
		words = line.split()
		for word in words:
			if word == '<ParameterRange': 
				labels.append(words[words.index(word)+1].split('"')[1])
				if words[words.index(word)+2].split('"')[1] == "length":Units.append("pc")
				if words[words.index(word)+2].split('"')[1] == "posangle":Units.append("deg")
				if words[words.index(word)+2].split('"')[1] == "mass":Units.append("Msun")		
				if words[words.index(word)+2].split('"')[1] == "dimless":Units.append("")				
fski.close()
print "\nFound following labels:"
for label in labels: print label

# fetching wavelengths from skifile
ski = open(skifile,'r')
nstellarcomps=0
for line in ski:
	if '<OligoStellarComp' in line: nstellarcomps+=1;	
	if '<OligoWavelengthGrid wavelengths' in line:
		wavelengths = (line.split('"')[1]+",").split(" micron,")
		wavelengths.pop()
print "\nNumber of stellar components: "+str(nstellarcomps)

# scan for best_simulation files
for root, dirnames, filenames in os.walk('.'):
	if glob.glob('./'+fskiname+'_BESTsimulations.dat') != []:
		matches.append(glob.glob('./'+fskiname+'_BESTsimulations.dat'))
	for dir in dirnames:
		if glob.glob("./"+dir+'/'+fskiname+'_BESTsimulations.dat') != []:
			matches.append(glob.glob("./"+dir+'/'+fskiname+'_BESTsimulations.dat'))
# reading last line of best_simulations while checking for inconsistencies
# before writing the solution file for each run
sol_length = 0
lowest_chi_simul = "" 
lowest_chi = 1e15
for match in matches:
	for bestsimulation in match:
		simulation = open(bestsimulation, 'r')
		for line in simulation:
			pass
		split_line = line.split(" ")
		split_line.pop()
		solutions.append(split_line)
		if sol_length == 0:	sol_length = len(split_line)			
		if sol_length != len(split_line):
			raise Exception("Inconsistency in number of free parameters between best_simulations: " +bestsimulation)
		if sol_length != len(wavelengths*(nstellarcomps+1))+2+len(labels):
			raise Exception("Inconsistency in number of free parameters for: " +bestsimulation)
		if float(split_line[1+len(labels)]) < float(lowest_chi):
			lowest_chi = split_line[1+len(labels)]
			lowest_chi_simul = len(solutions)-1
# adding Luminosities to labels
best_solution=solutions[lowest_chi_simul]
best_solution.pop(0)
best_chi=best_solution[len(labels)]
best_solution.pop(len(labels))
for i in range(0,nstellarcomps):
	for wavelength in wavelengths:
		labels.append("Lum"+str(i)+" at "+str(wavelength))
		Units.append("")

# creating separate solution files
std2=[0]*len(labels)
print "\nCreating following result files:"
for solution in solutions:
	chi=lowest_chi
	if solutions.index(solution) != lowest_chi_simul:
		chi=solution[len(labels)-nstellarcomps*len(wavelengths)+1]
		solution.pop(0)
		solution.pop(len(labels)-nstellarcomps*len(wavelengths))
	print (matches[solutions.index(solution)][0]).split(".dat")[0]+"_result.dat"
	output = open((matches[solutions.index(solution)][0]).split(".dat")[0]+"_result.dat",'w')
	output.write("{0:20} {1:12}".format("Parameter","Best"))
	output.write("\n------------------------------------------------------------------\n\n")
	for label in labels:
		output.write("{0:20} {1:12}".format(label, solution[labels.index(label)]))
		output.write("\n")
		if solutions.index(solution) != lowest_chi_simul:		
			std2[labels.index(label)]+= 1.0/(float(len(solutions)-1))*\
			(float(solution[labels.index(label)])-float(best_solution[labels.index(label)]))**2
	output.write("\nLowest chi: "+chi+"\n\n------------------------------------------------------------------")
	output.close()

# calculating and creating combined results
for value in std2:std2[std2.index(value)]=math.sqrt(float(value))
print "\nCombined result:"
print fskiname+"_RESULT.dat\n"
output_best = open(fskiname+"_RESULT.dat",'w')
output_best.write("Best simulation was simulation: "+matches[lowest_chi_simul][0]+"\n \n")
output_best.write("{0:20} {1:12} {2:3} {3:12} {4:8}".format("Parameter","Best","+-","Std","units"))
output_best.write("\n------------------------------------------------------------------------\n\n")
for label in labels:
	output_best.write("{0:20} {1:12} {2:3} {3:12} {4:8}".format(label,
	round(float(best_solution[labels.index(label)]),2), "+-", round(std2[labels.index(label)],2),Units[labels.index(label)]))
	output_best.write("\n")
output_best.write("\nLowest chi: "+lowest_chi+"\n\n------------------------------------------------------------------")
output_best.close()

# ------------------------------------------------------------------
