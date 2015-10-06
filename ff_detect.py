#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.3"
parser = argparse.ArgumentParser(prog='ff_detect', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@gmail.com)
git: https://github.com/jhelie/ff_detect
**********************************************
	
[ DESCRITPION ]

This script identifies the phospholipids which flip-flopped between two frames and
outputs a file with their respective selection strings following the format:
 -> 'resname,resid,starting_leaflet,beadname'
	

[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib


[ NOTES ]

1. It's a good idea to trjconv the xtc first and only outputs the relevant bead of each
   phospholipid as the script will run MUCH faster.	

2. Identification of the bilayer leaflets is controlled via 3 options:
   (a) selection of particles
    By default, the PO4 group is used to detect lipids and assess their flip-flop status.
    This bead name can be changed via the --bead flag, the selection string being:
    -> "name " + beadname
   
    Note that only lipids which contain the bead mentioned in the selection string
    will be taken into account to identify leaflets.
        
   (b) leaflet finding method: reference file
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use the default 15 Angstrom cutoff
    directly (without optimising):
     -> '--leaflets 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflets large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the gro file supplied in
    order to get a meaningful outcome.

3. Identification of lipids having flip-flopped can happen in two ways:
   (a) comparison of the results of the leaflets identification process
    By default leaflets are identified in the second structure file using the same method as
    in the first. Flip-flopping lipids are simply detected by comparing the content of the
    upper and lower leaflets in each file.

   (b) for systems in which it is not easy to identify leaflets in the second structure file
   (eg very large sytems deforming significantly or systems in which lipids in the process
   of flip-flopping 'bridge' the two leaflets together), the --neighbours option must be used.
   
   In this case the neighbouring lipds of each lipid are calculated and if more than half the
   neighbours of a given lipid belong to the opposite leaflet than its initial leaflet then
   this lipid is considered to have flip-flopped.
   
   The argument of --neighbours correspond to the max distance (in Angstrom) within which to
   consider neighbours. This method can be VERY slow for VERY large systems.


[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: reference structure [.gro] (required)
-g			: structure in which to identify flip-flops [.gro]
-o			: name of output folder

Lipids identification (see note 2)
-----------------------------------------------------
--bead		[PO4]	: lipids bead name, see note 2(a)
--leaflets		: leaflet identification, see note 2(b)
--neighbours	[15]	: flip-flopping lipids detection method, see note 3(b)
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='reffilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-g', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)

#lipids identification
parser.add_argument('--bead', nargs=1, dest='beadname', default=['PO4'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)
parser.add_argument('--neighbours', nargs='?', dest='neighbours', default="no", const="15", help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

#parse user inputs
#-----------------
args = parser.parse_args()
#data options
args.reffilename = args.reffilename[0]
args.grofilename = args.grofilename[0]
args.output_folder = args.output_folder[0]
#lipids identification
args.beadname = args.beadname[0]
args.cutoff_leaflet = args.cutoff_leaflet[0]

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.selections
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.reffilename):
	print "Error: file " + str(args.reffilename) + " not found."
	sys.exit(1)
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.neighbours != "no":
	args.neighbours = float(args.neighbours)
	if args.neighbours < 0:
		print "Error: --neighbours should be greater than 0, see note 3(b)"
		sys.exit(1)
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should be a number or 'large', see note 2"
		sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder=="no":
	args.output_folder="ff_detect_" + args.reffilename[:-4] + '_' + args.grofilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	#create folders
	os.mkdir(args.output_folder)

	#create log
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/ff_detect.log'
	output_log=open(filename_log, 'w')		
	output_log.write("[ff_detect v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python ff_detect.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

global nb_ff, nb_u2l, nb_l2u
global U_ref, U_gro
global upper_ref, upper_gro, lower_ref, lower_gro
nb_ff = 0
nb_u2l = 0
nb_l2u = 0
upper_to_lower = {}
lower_to_upper = {}

def data_loading():
	
	global U_ref, U_gro
	
	print "\nLoading files..."
	print " -" + str(args.reffilename) + "..."
	U_ref = Universe(args.reffilename)
	print " -" + str(args.grofilename) + "..."
	U_gro = Universe(args.grofilename)
	return

#=========================================================================================
# core functions
#=========================================================================================

def identify_leaflets_ref():

	global upper_ref
	global lower_ref
	
	print "\nProcessing reference file..."	
	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			tmp_cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U_ref, "name " + str(args.beadname))
			print " -identifying leaflets..."
			L_ref = MDAnalysis.analysis.leaflet.LeafletFinder(U_ref, "name " + str(args.beadname), tmp_cutoff_value[0])
		else:
			print " -identifying leaflets..."
			L_ref = MDAnalysis.analysis.leaflet.LeafletFinder(U_ref, "name " + str(args.beadname), args.cutoff_leaflet)
		if np.shape(L_ref.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		else:
			if L_ref.group(0).centerOfGeometry()[2] > L_ref.group(1).centerOfGeometry()[2]:
				upper_ref = L_ref.group(0)
				lower_ref = L_ref.group(1)
			else:
				upper_ref = L_ref.group(1)
				lower_ref = L_ref.group(0)

	#use cog 
	else:
		print " -identifying leaflets..."
		tmp_leaflets = U_ref.selectAtoms("name " + str(args.beadname))
		tmp_lipids_avg_z = tmp_leaflets.centerOfGeometry()[2]
		upper_ref = tmp_leaflets.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		lower_ref = tmp_leaflets.selectAtoms("prop z < " + str(tmp_lipids_avg_z))
	
	#display results	
	print " -found 2 leaflets: ", upper_ref.numberOfResidues(), "(upper) and ", lower_ref.numberOfResidues(), "(lower) lipids"

	return
def identify_leaflets_gro():
	
	global upper_gro
	global lower_gro

	print "\nProcessing second file..."
	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			tmp_cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U_gro, "name " + str(args.beadname))
			print " -identifying leaflets..."
			L_gro = MDAnalysis.analysis.leaflet.LeafletFinder(U_gro, "name " + str(args.beadname), tmp_cutoff_value[0])
		else:
			print " -identifying leaflets..."
			L_gro = MDAnalysis.analysis.leaflet.LeafletFinder(U_gro, "name " + str(args.beadname), args.cutoff_leaflet)
		if np.shape(L_gro.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		else:
			if L_gro.group(0).centerOfGeometry()[2] > L_gro.group(1).centerOfGeometry()[2]:
				upper_gro = L_gro.group(0)
				lower_gro = L_gro.group(1)
			else:
				upper_gro = L_gro.group(1)
				lower_gro = L_gro.group(0)

	#use cog 
	else:
		print " -identifying leaflets..."
		tmp_leaflets = U_gro.selectAtoms("name " + str(args.beadname))
		tmp_lipids_avg_z = tmp_leaflets.centerOfGeometry()[2]
		upper_gro = tmp_leaflets.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		lower_gro = tmp_leaflets.selectAtoms("prop z < " + str(tmp_lipids_avg_z))


	#display results	
	print " -found 2 leaflets: " + str(upper_gro.numberOfResidues()) + "(upper) and " + str(lower_gro.numberOfResidues()) + "(lower) lipids"
	
	return
def identify_ff():

	global nb_ff, nb_u2l, nb_l2u
	
	print "\nIdentifying flip-flopping lipids..."
	
	#method 1: compare leaflets
	#--------------------------
	if args.neighbours == "no":
		
		#upper to lower: lipids in upper_ref NOT in upper_gro
		upref_to_lwgro = np.setdiff1d(upper_ref.resnums(),upper_gro.resnums())
		for r_num in upref_to_lwgro:
			tmp_a = upper_ref.selectAtoms("resid " + str(r_num))
			upper_to_lower[tmp_a.resnums()[0]] = tmp_a.resnames()[0]
		
		#lower to upper: lipids in lower_ref NOT in lower_gro
		lwref_to_upgro = np.setdiff1d(lower_ref.resnums(),lower_gro.resnums())
		for r_num in lwref_to_upgro:
			tmp_a = lower_ref.selectAtoms("resid " + str(r_num))
			lower_to_upper[tmp_a.resnums()[0]] = tmp_a.resnames()[0]
		
	#method 2: compare neighbours
	#----------------------------
	else:
		tmp_upper_neighbours = {}
		tmp_lower_neighbours = {}
		tmp_upper_nb = upper_ref.numberOfResidues()
		tmp_lower_nb = lower_ref.numberOfResidues()
		tmp_upper_resnums = upper_ref.resnums()
		tmp_lower_resnums = lower_ref.resnums()
		tmp_U_gro_lip = U_gro.selectAtoms("name " + str(args.beadname))
		
		#upper to lower
		for a_index in range(0,tmp_upper_nb):
			#display progress
			progress = '\r -upper leaflet: processing lipid ' + str(a_index+1) + '/' + str(tmp_upper_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			
			#check whether more than half of the neighbours belong to the lower leaflet
			tmp_a = upper_ref[a_index]
			tmp_neighbours = tmp_U_gro_lip.selectAtoms("around " + str(args.neighbours) + " resid " + str(tmp_a.resid)).resnums()
			tmp_neighbours = np.in1d(tmp_neighbours, tmp_lower_resnums)
			if len(tmp_neighbours) == 0:
				print ""
				print "  Warning: no lipid neighbours within", str(args.neighbours), "of bead", str(args.beadname), "of", str(tmp_a.resname), str(tmp_a.resid) 
			else:
				tmp_upper_neighbours[a_index] =  len(tmp_neighbours[tmp_neighbours==True]) / len(tmp_neighbours)
				if tmp_upper_neighbours[a_index] > 0.5:
					upper_to_lower[tmp_a.resnum] = tmp_a.resname
		print ''
	
		#lower to upper
		for a_index in range(0,tmp_lower_nb):
			#display progress
			progress = '\r -lower leaflet: processing lipid ' + str(a_index+1) + '/' + str(tmp_lower_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			
			#check whether more than half of the neighbours belong to the lower leaflet
			tmp_a = lower_ref[a_index]
			tmp_neighbours = tmp_U_gro_lip.selectAtoms("around " + str(args.neighbours) + " resid " + str(tmp_a.resid)).resnums()
			tmp_neighbours = np.in1d(tmp_neighbours, tmp_upper_resnums)
			if len(tmp_neighbours) == 0:
				print ""
				print "  Warning: no lipid neighbours within ", str(args.neighbours), "of bead ", str(args.beadname), "of ", str(tmp_a.resname), str(tmp_a.resid) 
				print ""
			else:
				tmp_lower_neighbours[a_index] =  len(tmp_neighbours[tmp_neighbours==True]) / float(len(tmp_neighbours))
				if tmp_lower_neighbours[a_index] > 0.5:				
					lower_to_upper[tmp_a.resnum] = tmp_a.resname
		print ''
		
	#count nb of flip-flopping lipids
	#--------------------------------
	nb_u2l = len(upper_to_lower.keys())
	nb_l2u = len(lower_to_upper.keys())
	nb_ff = nb_u2l + nb_l2u
		
	return

#=========================================================================================
# outputs
#=========================================================================================

def write_sele_file():
	
	#open file
	print " - .sele file..."
	output_sele = open(os.getcwd() + '/' + args.output_folder + '/flipflop.sele', 'w')

	#upper to lower
	#--------------
	if nb_u2l > 0:
		tmp_sorted = sorted(upper_to_lower.items(), key=operator.itemgetter(1))
		for r in range(0,nb_u2l):
			output_sele.write(str(tmp_sorted[r][1]) + "," + str(tmp_sorted[r][0]) + ",upper," + str(args.beadname) + "\n")
	
	#lower to upper
	#--------------
	if nb_l2u > 0:
		tmp_sorted = sorted(lower_to_upper.items(), key=operator.itemgetter(1))
		for r in range(0,nb_l2u):
			output_sele.write(str(tmp_sorted[r][1]) + "," + str(tmp_sorted[r][0]) + ",lower," + str(args.beadname) + "\n")
	
	output_sele.close()
	
	return
def write_pml_vmd_file():

	print " - .pml and .vmd files..."
	
	#upper to lower
	#--------------
	if nb_u2l > 0:
		#create writers instances
		leaflet_writer_vmd = MDAnalysis.selections.vmd.SelectionWriter(args.output_folder + "/upper")
		leaflet_writer_pml = MDAnalysis.selections.pymol.SelectionWriter(args.output_folder + "/upper")
		ff_writer_vmd = MDAnalysis.selections.vmd.SelectionWriter(args.output_folder + "/ff_upper_to_lower")
		ff_writer_pml = MDAnalysis.selections.pymol.SelectionWriter(args.output_folder + "/ff_upper_to_lower")

		#create selection string
		for r_index in range(0,len(upper_to_lower.keys())):
			r_num = upper_to_lower.keys()[r_index]
			if r_index == 0:
				tmp_ff_u2l_string = "(resname " + str(upper_to_lower[r_num]) + " and resnum " + str(r_num) + ")"
			else:
				tmp_ff_u2l_string += " or (resname " + str(upper_to_lower[r_num]) + " and resnum " + str(r_num) + ")"
		
		#create selections
		upper_clean = upper_ref.selectAtoms(" not (" + str(tmp_ff_u2l_string) + ")")
		tmp_u2l = U_ref.selectAtoms(str(tmp_ff_u2l_string))
		
		#write selections
		leaflet_writer_vmd.write(upper_clean, name = "upper")
		leaflet_writer_pml.write(upper_clean, name = "upper")
		ff_writer_vmd.write(tmp_u2l, name = "ff_u2l")
		ff_writer_pml.write(tmp_u2l, name = "ff_u2l")
	else:
		leaflet_writer_vmd = MDAnalysis.selections.vmd.SelectionWriter(args.output_folder + "/upper")
		leaflet_writer_pml = MDAnalysis.selections.pymol.SelectionWriter(args.output_folder + "/upper")
		leaflet_writer_vmd.write(upper_ref, name = "upper")
		leaflet_writer_pml.write(upper_ref, name = "upper")
			
	#lower to upper
	#--------------
	if nb_l2u > 0:
		#create writers instances
		leaflet_writer_vmd = MDAnalysis.selections.vmd.SelectionWriter(args.output_folder + "/lower")
		leaflet_writer_pml = MDAnalysis.selections.pymol.SelectionWriter(args.output_folder + "/lower")
		ff_writer_vmd = MDAnalysis.selections.vmd.SelectionWriter(args.output_folder + "/ff_lower_to_upper")
		ff_writer_pml = MDAnalysis.selections.pymol.SelectionWriter(args.output_folder + "/ff_lower_to_upper")
		
		#create selection string
		for r_index in range(0,len(lower_to_upper.keys())):
			r_num = lower_to_upper.keys()[r_index]
			if r_index == 0:
				tmp_ff_l2u_string = "(resname " + str(lower_to_upper[r_num]) + " and resnum " + str(r_num) + ")"
			else:
				tmp_ff_l2u_string += " or (resname " + str(lower_to_upper[r_num]) + " and resnum " + str(r_num) + ")"
		
		#create selections
		lower_clean = lower_ref.selectAtoms(" not (" + str(tmp_ff_l2u_string) + ")")
		tmp_l2u = U_ref.selectAtoms(str(tmp_ff_l2u_string))
		
		#write selections
		leaflet_writer_vmd.write(lower_clean, name = "lower")
		leaflet_writer_pml.write(lower_clean, name = "lower")
		ff_writer_vmd.write(tmp_l2u, name = "ff_l2u")
		ff_writer_pml.write(tmp_l2u, name = "ff_l2u")
	else:
		leaflet_writer_vmd = MDAnalysis.selections.vmd.SelectionWriter(args.output_folder + "/lower")
		leaflet_writer_pml = MDAnalysis.selections.pymol.SelectionWriter(args.output_folder + "/lower")
		leaflet_writer_vmd.write(lower_ref, name = "lower")
		leaflet_writer_pml.write(lower_ref, name = "lower")

	return


##########################################################################################
# ALGORITHM
##########################################################################################

#load data
data_loading()

#identify leaflets in file(s) supplied
identify_leaflets_ref()
if args.neighbours == "no":
	identify_leaflets_gro()

#identifiy flip-flops
identify_ff()

#write outputs
if nb_ff > 0:
	print "\nWriting selection files..."
	write_sele_file()
	write_pml_vmd_file()
	print "\nFinished successfully!", str(nb_ff), "flip-flops detected (" + str(nb_u2l) + " u->l and " + str(nb_l2u) + " l->u), check results in folder '" + str(args.output_folder) + "'"
else:
	print "\nFinished successfully! 0 flip-flops detected."

#exit
sys.exit(0)
