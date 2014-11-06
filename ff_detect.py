# this script is called by flipflop.sh

#################################################################
# Import python modules
#################################################################
#import general python tools
import operator
from operator import itemgetter
import sys, os, shutil
import math

#import python extensions/packages to manipulate arrays
import numpy 				#to manipulate arrays
import scipy 				#mathematical tools and recipesimport MDAnalysis

#import graph building module
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
import matplotlib.cm				#colours library
import matplotlib.ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
fontP=FontProperties()

#import MDAnalysis
import MDAnalysis
from MDAnalysis import *
import MDAnalysis.analysis
import MDAnalysis.analysis.leaflet
import MDAnalysis.analysis.distances

#set MDAnalysis to use periodic boundary conditions
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

################################################################################################################################################
# RETRIEVE USER INPUTS
################################################################################################################################################

#ff_detect.py $version $reffilename $grofilename $xtcfilename $optimise_opt $outputname
version=sys.argv[1]
reffilename=sys.argv[2]
grofilename=sys.argv[3]
xtcfilename=sys.argv[4]
optimise_opt=sys.argv[5]
outputname=sys.argv[6]

################################################################################################################################################
# DATA LOADING
################################################################################################################################################

#data declatation
global other1_pres, other2_pres, other1_dict, other2_dict
global upper1, upper2, lower1, lower2
global ff_u1_to_l2, ff_l1_to_u2, ff_u2_to_l1, ff_l2_to_u1
global up1_to_lw2, lw1_to_up2, up2_to_lw1, lw2_to_up1
other1_pres="no"
other2_pres="no"
other1_dict={}
other2_dict={}
ff_u1_to_l2={}
ff_l1_to_u2={}
ff_u2_to_l1={}
ff_l2_to_u1={}

#load universe
if xtcfilename=="no":
	print "\nLoading files..."
	U1=Universe(reffilename)
	U2=Universe(grofilename)
else:
	print "\nLoading trajectory..."
	U1=Universe(reffilename, xtcfilename)
	tmp_all=U1.selectAtoms("not resname W and not resname ION")
	ts_start=U1.trajectory[0]
	ts_start_time=str(int(ts_start.time/float(1000)))
	tmp_all.write(os.getcwd() + '/' + outputname + '/xtc_start_' + ts_start_time + 'ns.gro')
	ts_end=U1.trajectory[U1.trajectory.numframes-1]
	ts_end_time=str(int(ts_end.time/float(1000)))
	tmp_all.write(os.getcwd() + '/' + outputname + '/xtc_end_' + ts_end_time + 'ns.gro')
	U1=Universe(os.getcwd() + '/' + outputname + '/xtc_start_' + ts_start_time + 'ns.gro')
	U2=Universe(os.getcwd() + '/' + outputname + '/xtc_end_' + ts_end_time + 'ns.gro')

################################################################################################################################################
# FUNCTIONS
################################################################################################################################################

def process_leaflets():

	global up1_to_lw2, lw1_to_up2, up2_to_lw1, lw2_to_up1
	global upper1, upper2, lower1, lower2
	
	#identify leaflets
	#-----------------
	if optimise_opt=="yes":
		print " -optimising cutoff for file1..."
		cutoff_value1=MDAnalysis.analysis.leaflet.optimize_cutoff(U1, "name PO4 or name PO3 or name B1A") #name PO4 or name AM1 or name PO1

		#debug
		print cutoff_value1

		print " -identifying leaflets for file1..."
		L1=MDAnalysis.analysis.leaflet.LeafletFinder(U1, "name PO4 or name PO3 or name B1A", cutoff_value1[0])

		print " -optimising cutoff for file2..."
		cutoff_value2=MDAnalysis.analysis.leaflet.optimize_cutoff(U2, "name PO4 or name PO3 or name B1A")
		print " -identifying leaflets for file2..."
		L2=MDAnalysis.analysis.leaflet.LeafletFinder(U2, "name PO4 or name PO3 or name B1A", cutoff_value2[0])
	else:
		print " -identifying leaflets for file1..."
		L1=MDAnalysis.analysis.leaflet.LeafletFinder(U1, "name PO4 or name PO3 or name B1A", 11)
		print " -identifying leaflets for file2..."
		L2=MDAnalysis.analysis.leaflet.LeafletFinder(U2, "name PO4 or name PO3 or name B1A", 11)
		
	#debug
	print L1.groups()
	print L2.groups()

	#check leaflets
	#--------------
	#gro file 1
	if numpy.shape(L1.groups())[0]<2:
		print "Error: imposssible to identify 2 leaflets."
		sys.exit(1)
	else:
		if L1.group(0).centerOfGeometry()[2] > L1.group(1).centerOfGeometry()[2]:
			upper1=L1.group(0)
			lower1=L1.group(1)
		else:
			upper1=L1.group(1)
			lower1=L1.group(0)
	
		if numpy.shape(L1.groups())[0]>2:
			other1_pres="yes"
			other1=L1.group(2)
		
			if numpy.shape(L1.groups())[0]>3:
				for g in range(3,numpy.shape(L1.groups())[0]):
					other1+=L1.group(g)
					
			for r in other1.resnums():
				tmp=U1.selectAtoms("resid " + str(r))
				if tmp.resnames()[0] in other1_dict.keys():
					other1_dict[tmp.resnames()[0]].append(r)
				else:
					other1_dict[tmp.resnames()[0]]=[]
					other1_dict[tmp.resnames()[0]].append(r)
					
	#gro file 2
	if numpy.shape(L2.groups())[0]<2:
		print "Error: imposssible to identify 2 leaflets."
		sys.exit(1)
	else:	
		if L2.group(0).centerOfGeometry()[2] > L2.group(1).centerOfGeometry()[2]:
			upper2=L2.group(0)
			lower2=L2.group(1)
		else:
			upper2=L2.group(1)
			lower2=L2.group(0)
	
		if numpy.shape(L2.groups())[0]>2:
			other2_pres="yes"
			other2=L2.group(2)
		
			if numpy.shape(L2.groups())[0]>3:
				for g in range(3,numpy.shape(L2.groups())[0]):
					other2+=L2.group(g)
					
			for r in other2.resnums():
				tmp=U1.selectAtoms("resid " + str(r))
				if tmp.resnames()[0] in other2_dict.keys():
					other2_dict[tmp.resnames()[0]].append(r)
				else:
					other2_dict[tmp.resnames()[0]]=[]
					other2_dict[tmp.resnames()[0]].append(r)
	
	#display
	#-------
	print " -file1: found " + str(upper1.numberOfResidues()) + " (upper) and " + str(lower1.numberOfResidues()) + " (lower) lipids"
	print " -file2: found " + str(upper2.numberOfResidues()) + " (upper) and " + str(lower2.numberOfResidues()) + " (lower) lipids"
	
	#compare leaflets
	#----------------
	#lipids in upper1 NOT in upper2
	up1_to_lw2=numpy.setdiff1d(upper1.resnums(),upper2.resnums())
	#lipids in lower1 NOT in lower2
	lw1_to_up2=numpy.setdiff1d(lower1.resnums(),lower2.resnums())
	#lipids in upper2 NOT in upper1
	up2_to_lw1=numpy.setdiff1d(upper2.resnums(),upper1.resnums())
	#lipids in lower2 NOT in lower1
	lw2_to_up1=numpy.setdiff1d(lower2.resnums(),lower1.resnums())
	
	return
def calc_statistics():
	
	global ff_u1_to_l2, ff_l1_to_u2, ff_u2_to_l1, ff_l2_to_u1
	global up1_to_lw2, lw1_to_up2, up2_to_lw1, lw2_to_up1

	#lipids in upper1 NOT in upper2
	for r in up1_to_lw2:
		tmp=U1.selectAtoms("resid " + str(r))
		if tmp.resnames()[0] in ff_u1_to_l2.keys():
			ff_u1_to_l2[tmp.resnames()[0]].append(r)
		else:
			ff_u1_to_l2[tmp.resnames()[0]]=[]
			ff_u1_to_l2[tmp.resnames()[0]].append(r)
	
	#lipids in lower1 NOT in lower2
	for r in lw1_to_up2:
		tmp=U1.selectAtoms("resid " + str(r))
		if tmp.resnames()[0] in ff_l1_to_u2.keys():
			ff_l1_to_u2[tmp.resnames()[0]].append(r)
		else:
			ff_l1_to_u2[tmp.resnames()[0]]=[]
			ff_l1_to_u2[tmp.resnames()[0]].append(r)
	
	#lipids in upper2 NOT in upper1
	for r in up2_to_lw1:
		tmp=U1.selectAtoms("resid " + str(r))
		if tmp.resnames()[0] in ff_u2_to_l1.keys():
			ff_u2_to_l1[tmp.resnames()[0]].append(r)
		else:
			ff_u2_to_l1[tmp.resnames()[0]]=[]
			ff_u2_to_l1[tmp.resnames()[0]].append(r)
	
	#lipids in lower2 NOT in lower1
	for r in lw2_to_up1:
		tmp=U1.selectAtoms("resid " + str(r))
		if tmp.resnames()[0] in ff_l2_to_u1.keys():
			ff_l2_to_u1[tmp.resnames()[0]].append(r)
		else:
			ff_l2_to_u1[tmp.resnames()[0]]=[]
			ff_l2_to_u1[tmp.resnames()[0]].append(r)

	return
def write_ff_stats():

	global ff_u1_to_l2, ff_l1_to_u2, ff_u2_to_l1, ff_l2_to_u1
	global upper1, upper2, lower1, lower2
	
	output_stat = open(os.getcwd() + '/' + outputname + '/flipflop.stat', 'w')
	output_stat.write("This file was generated by the script flipflop v" + str(version) +"\n")
	output_stat.write("\n")
	output_stat.write("Leaflets identification\n")
	output_stat.write("=======================\n")
	output_stat.write("\n")
	output_stat.write("      	File 1		File 2\n")
	output_stat.write("      	------		------\n")
	output_stat.write("upper:	" + str(upper1.numberOfResidues()) + "		" + str(upper2.numberOfResidues()) + "\n")
	output_stat.write("lower:	" + str(lower1.numberOfResidues()) + "		" + str(lower2.numberOfResidues()) + "\n")
	if other1_pres=="yes" and other2_pres=="yes":
		output_stat.write("other:	" + str(other1.numberOfResidues()) + "		" + str(other2.numberOfResidues()) + "\n")
		output_stat.write("\n")
		output_stat.write("other file 1:\n")
		for r in other1_dict.keys():
			output_stat.write("-" + str(r) + ":" + str(numpy.size(other1_dict[r])) + " " + str(other1_dict[r]) + "\n")
		output_stat.write("\n")
		output_stat.write("other file 2:\n")
		for r in other2_dict.keys():
			output_stat.write("-" + str(r) + ":" + str(numpy.size(other2_dict[r])) + " " + str(other2_dict[r]) + "\n")
	elif other1_pres=="yes" and other2_pres=="no":
		output_stat.write("other:	" + str(other1.numberOfResidues()) + "\n")
		output_stat.write("\n")
		output_stat.write("other file 1:\n")
		for r in other1_dict.keys():
			output_stat.write("-" + str(r) + ":" + str(numpy.size(other1_dict[r])) + " " + str(other1_dict[r]) + "\n")
	elif other1_pres=="no" and other2_pres=="yes":
		output_stat.write("other:	      		" + str(other2.numberOfResidues()) + "\n")
		output_stat.write("\n")
		output_stat.write("other file 2:\n")
		for r in other2_dict.keys():
			output_stat.write("-" + str(r) + ":" + str(numpy.size(other2_dict[r])) + " " + str(other2_dict[r]) + "\n")
	output_stat.write("\n")
	output_stat.write("\n")	
	output_stat.write("Flip-flops statistics\n")
	output_stat.write("=====================\n")
	output_stat.write("\n")	
	output_stat.write("upper 1 -> lower 2:\n")
	output_stat.write("-------------------\n")
	for r in ff_u1_to_l2.keys():
		output_stat.write("-" + str(r) + ":" + str(numpy.size(ff_u1_to_l2[r])) + " " + str(ff_u1_to_l2[r]) + "\n")
	output_stat.write("\n")	
	output_stat.write("lower 1 -> upper 2:\n")
	output_stat.write("-------------------\n")
	for r in ff_l1_to_u2.keys():
		output_stat.write("-" + str(r) + ":" + str(numpy.size(ff_l1_to_u2[r])) + " " + str(ff_l1_to_u2[r]) + "\n")
	output_stat.write("\n")	
	output_stat.write("\n")	
	output_stat.write("Consistency check\n")
	output_stat.write("=================\n")
	output_stat.write("\n")	
	output_stat.write("If the figures below and above are different there is probably some issue in identifying leaflets in one of the files provided.\n")	
	output_stat.write("\n")	
	output_stat.write("upper 2 -> lower 1:\n")
	output_stat.write("-------------------\n")
	for r in ff_u2_to_l1.keys():
		output_stat.write("-" + str(r) + ":" + str(numpy.size(ff_u2_to_l1[r])) + " " + str(ff_u2_to_l1[r]) + "\n")
	output_stat.write("\n")	
	output_stat.write("lower 2 -> upper 1:\n")
	output_stat.write("-------------------\n")
	for r in ff_l2_to_u1.keys():
		output_stat.write("-" + str(r) + ":" + str(numpy.size(ff_l2_to_u1[r])) + " " + str(ff_l2_to_u1[r]) + "\n")
	output_stat.write("\n")	
	output_stat.close()

	return
def write_selection_file():
	#-------------------------------------------------------------------
	# Write selection file which can be used as input by flipflop_stat
	#-------------------------------------------------------------------
	global ff_u1_to_l2, ff_l1_to_u2, ff_u2_to_l1, ff_l2_to_u1	
	
	output_sele = open(os.getcwd() + '/' + outputname + '/flipflop.sele', 'w')
	for s in ff_u1_to_l2:
		for r in ff_u1_to_l2[s]:
			output_sele.write(str(s) + "," + str(r) + ",upper,PO4\n")
	for s in ff_l1_to_u2:
		for r in ff_l1_to_u2[s]:
			output_sele.write(str(s) + "," + str(r) + ",lower,PO4\n")
	output_sele.close()
	
	return
def graph_ff_stats():

	#-------------------------------------------------------------------
	#-what: nb of flipflopping lipids
	#-plot: bar chart, 2 bars (outer2inner and inner2outer) for each specie
	#-------------------------------------------------------------------		
	global ff_u1_to_l2, ff_l1_to_u2, ff_u2_to_l1, ff_l2_to_u1
	
	filename_png=os.getcwd() + '/' + outputname + '/flipflops.png'
	filename_svg=os.getcwd() + '/' + outputname + '/flipflops.svg'
	
	#create figure
	#-------------
	fig=plt.figure(figsize=(3.25, 3)) 								#1 column format
	xticks_pos=numpy.arange(1,4)
	xticks_lab=['POPC', 'POPE', 'POPS']

	#create data
	#-----------
	#upper to lower
	if "POPC" in ff_u1_to_l2.keys():
		nb_POPC_upper2lower=numpy.size(ff_u1_to_l2["POPC"])
	else:
		nb_POPC_upper2lower=0
	if "POPE" in ff_u1_to_l2.keys():
		nb_POPE_upper2lower=numpy.size(ff_u1_to_l2["POPE"])
	else:
		nb_POPE_upper2lower=0
	if "POPS" in ff_u1_to_l2.keys():
		nb_POPS_upper2lower=numpy.size(ff_u1_to_l2["POPS"])
	else:
		nb_POPS_upper2lower=0
	#lower to upper
	if "POPC" in ff_l1_to_u2.keys():
		nb_POPC_lower2upper=numpy.size(ff_l1_to_u2["POPC"])
	else:
		nb_POPC_lower2upper=0
	if "POPE" in ff_l1_to_u2.keys():
		nb_POPE_lower2upper=numpy.size(ff_l1_to_u2["POPE"])
	else:
		nb_POPE_lower2upper=0
	if "POPS" in ff_l1_to_u2.keys():
		nb_POPS_lower2upper=numpy.size(ff_l1_to_u2["POPS"])
	else:
		nb_POPS_lower2upper=0
		
	#plot data: nb of flipflops
	#-------------------------
	ax1 = fig.add_subplot(111)
	plt.bar(xticks_pos-0.25, [nb_POPC_upper2lower, nb_POPE_upper2lower, nb_POPS_upper2lower], width=0.25, color='k', label='upper -> lower')
	plt.bar(xticks_pos-0.00, [nb_POPC_lower2upper, nb_POPE_lower2upper, nb_POPS_lower2upper], width=0.25, color='w', label='lower -> upper', hatch='/')
	plt.xticks(xticks_pos, xticks_lab)
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("nb of flipflops", size='medium')

	#legend, labels, etc
	#-------------------
	#ax1.set_ylabel('nb of peptides', fontsize="small")
	#ax2.set_ylabel('% of peptides', fontsize="small")
	ax1.tick_params(axis='both', direction='out')
	ax1.spines["right"].set_visible(False)								# remove unneeded axes
	ax1.spines["top"].set_visible(False)
	ax1.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax1.get_yaxis().tick_left()
	ax1.set_xlim(0.5, 3.5)
	ax1.set_ylim(ymin=0)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small")
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax1.yaxis.get_major_formatter().set_powerlimits((-3, 3))

	#save figure
	#-----------
	plt.subplots_adjust(left=0.1, right=0.96)
	fig.savefig(filename_svg)
	fig.savefig(filename_png)
	plt.close()

	return

################################################################################################################################################
# ALGORITHM
################################################################################################################################################

#identify leaflets
print "\nIdentifying leaflets..."
process_leaflets()

#identify flip-flopping lipids
print "\nIdentifying flip-flopping lipids..."
calc_statistics()
nb_ff=numpy.size(up1_to_lw2) + numpy.size(lw1_to_up2)
print " -found " + str(nb_ff) + " flipflopping lipids."

#write outputs
if nb_ff>0:
	#on screen display
	u2l_results=""
	l2u_results=""
	for r in ff_u1_to_l2.keys():
		u2l_results+=" " + str(r) + " (" + str(numpy.size(ff_u1_to_l2[r])) + ")"
	for r in ff_l1_to_u2.keys():
		l2u_results+=" " + str(r) + " (" + str(numpy.size(ff_l1_to_u2[r])) + ")"
	print " -upper to lower: "  + u2l_results
	print " -lower to upper: "  + l2u_results
	
	#write files
	print "\nWriting outputs..."
	write_ff_stats()
	write_selection_file()
	graph_ff_stats()

#exit
sys.exit(0)
