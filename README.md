ff_detect
=========

Python script to detect phospholipids flip-flops between two frames.
To print the help below: ```python ff_detect.py --help```

```
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
```
