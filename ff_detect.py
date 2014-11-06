#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

################################################################################################################################################
# RETRIEVE USER INPUTS
################################################################################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog='ff_detect', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
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

   (b) for systems in which lipids are still flip-flopping and the bilayer deforms significantly
   an other approach can be specified via the --neighbours flag.
   In this case the neighbouring lipds of each lipid are calculated and if more than half the
   neighbours of a given lipid belong to the opposite leaflet than its initial leaflet then this
   lipid is considered to have flip-flopped.
   For large systems deforming a lot and involving several flip-flpos it is, for now, the only
   solution but it is VERY slow. But you should only have to do it once to get the list of
   flip-flopping lipids.
   The argument of --neighbours correspond to the max distance (in Angstrom) within which to
   consider neighbours.

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
parser.add_argument('--neighbours', nargs='?', dest='neighbours', default=["no"], const=[15], help=argparse.SUPPRESS)

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
args.neighbours = args.neighbours[0]

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
	args.output_folder="ff_" + args.reffilename[:-4] + '_' + args.grofilename[:-4]
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

global U_ref, U_gro
global upper_ref, upper_gro, lower_ref, lower_gro

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


def data_loading():
	
	global U_ref, U_gro
	
	print "\nLoading files..."
	print " -file1..."
	U_ref = Universe(reffilename)
	print " -file2..."
	U_gro = Universe(grofilename)
	return
def identify_leaflets_ref():

	global upper_ref
	global lower_ref
	
	print "\nProcessing reference file..."	
	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff for file1..."
			tmp_cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U_ref, "name " + str(args.beadname))
			print " -identifying leaflets for file1..."
			L_ref = MDAnalysis.analysis.leaflet.LeafletFinder(U_ref, "name " + str(args.beadname), tmp_cutoff_value[0])
		else:
			print " -identifying leaflets for file1..."
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
		print " -identifying leaflets for file1..."
		tmp_leaflets = U_ref.selectAtoms("name " + str(args.beadname))
		tmp_lipids_avg_z = tmp_leaflets.centerOfGeometry()[2]
		upper_ref = tmp_leaflets.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		lower_ref = tmp_leaflets.selectAtoms("prop z < " + str(tmp_lipids_avg_z))
	
	#display results	
	print " -found 2 leaflets: " + upper_ref.numberOfResidues() + "(upper) and " + lower_ref.numberOfResidues() + "(lower) lipids"

	return
def identify_leaflets_gro():
	
	global upper_gro
	global lower_gro

	print "\nProcessing second file..."
	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff for file1..."
			tmp_cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U_gro, "name " + str(args.beadname))
			print " -identifying leaflets for file1..."
			L_gro = MDAnalysis.analysis.leaflet.LeafletFinder(U_gro, "name " + str(args.beadname), tmp_cutoff_value[0])
		else:
			print " -identifying leaflets for file1..."
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
		print " -identifying leaflets for file1..."
		tmp_leaflets = U_gro.selectAtoms("name " + str(args.beadname))
		tmp_lipids_avg_z = tmp_leaflets.centerOfGeometry()[2]
		upper_gro = tmp_leaflets.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		lower_gro = tmp_leaflets.selectAtoms("prop z < " + str(tmp_lipids_avg_z))


	#display results	
	print " -found 2 leaflets: " + upper_gro.numberOfResidues() + "(upper) and " + lower_gro.numberOfResidues() + "(lower) lipids"
	
	return

def identify_ff():

	#method 1: compare leaflets
	#--------------------------
	if args.neighbours == "no":
		
		#compare contents
		upref_to_lwgro = np.setdiff1d(upper_ref.resnums(),upper_gro.resnums())		#lipids in upper_ref NOT in upper_gro
		lwref_to_upgro = np.setdiff1d(lower_ref.resnums(),lower_gro.resnums())		#lipids in lower_ref NOT in lower_gro
		upgro_to_lwref = np.setdiff1d(upper_gro.resnums(),upper_ref.resnums())		#lipids in upper_gro NOT in upper_ref
		lwgro_to_upref = np.setdiff1d(lower_gro.resnums(),lower_ref.resnums())		#lipids in lower_gro NOT in lower_ref

		#build lists
		for r in upref_to_lwgro:										#lipids in upper1 NOT in upper2
			tmp = U_ref.selectAtoms("resid " + str(r))
			if tmp.resnames()[0] in ff_uref_to_lgro.keys():
				ff_uref_to_lgro[tmp.resnames()[0]].append(r)
			else:
				ff_uref_to_lgro[tmp.resnames()[0]]=[]
				ff_uref_to_lgro[tmp.resnames()[0]].append(r)
		for r in lwref_to_upgro:										#lipids in lower1 NOT in lower2
			tmp = U_ref.selectAtoms("resid " + str(r))
			if tmp.resnames()[0] in ff_lref_to_ugro.keys():
				ff_lref_to_ugro[tmp.resnames()[0]].append(r)
			else:
				ff_lref_to_ugro[tmp.resnames()[0]]=[]
				ff_lref_to_ugro[tmp.resnames()[0]].append(r)
		for r in upgro_to_lwref:										#lipids in upper2 NOT in upper1
			tmp=U_ref.selectAtoms("resid " + str(r))
			if tmp.resnames()[0] in ff_ugro_to_lref.keys():
				ff_ugro_to_lref[tmp.resnames()[0]].append(r)
			else:
				ff_ugro_to_lref[tmp.resnames()[0]]=[]
				ff_ugro_to_lref[tmp.resnames()[0]].append(r)
		for r in lwgro_to_upref:										#lipids in lower2 NOT in lower1
			tmp=U_ref.selectAtoms("resid " + str(r))
			if tmp.resnames()[0] in ff_lgro_to_uref.keys():
				ff_lgro_to_uref[tmp.resnames()[0]].append(r)
			else:
				ff_lgro_to_uref[tmp.resnames()[0]]=[]
				ff_lgro_to_uref[tmp.resnames()[0]].append(r)	

	#method 2: compare neighbours
	#----------------------------
	else:
		print "to do"

	
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
