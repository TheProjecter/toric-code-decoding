--------------------------------------------------------
      CONCATENATED [[4,2,2]]-TORIC CODE SIMULATOR

	Martin Suchara, Tobias Hoelzer, Ben Criger
	2012/September 2013/April 2014
	UC Berkeley/RWTH Aachen University
	
	See Bachelor Thesis: 
		Study of [[4,2,2]]-concatenated toric code
		for a scalable circuit-QED architecture
--------------------------------------------------------

VERSION DESCRIPTION:

Aids in the derivation of threshold estimates for the toric code vs. 
a depolarizing error channel, and a concatenated code defined on a
square-octagon lattice vs. the same model. 

"This version is only slightly different than the one from Martin Suchara.
Only Z and Y errors are removed." (Note: Find out what this means. --Ben)

1.1 COMPILATION

To compile type make clean; make.
To run the simulation type ./qecsim.

1.2 SIMULATION PARAMETERS

The simulation parameters are specified in in/parameters.txt. The file has the following format:

XSize = 2 4 6
iterations = 100000
pMin = 0.001
pMax = 0.1
pStep = 0.001
is_concatenated = 1
p_sqr_factor = 1.0
p_oct_factor = 1.0

These variables are:
 + *XSize* a space-separated integer list of lattice lengths to be 
   tested. For each N in XSize, an N-by-N lattice will be studied.
 + *iterations* an integer denoting how many independent trials 
   are to be carried out for each value of p (the physical error
   probability), for each lattice size. 
 + *pMin, pMax, pStep* define the linear space of physical error 
   probabilities to be examined for each lattice size. 
 + *is_concatenated* 0/1 value which indicates whether to use
   the concatenated code (1) or the vanilla toric code(0).
 + *p_sqr_factor, p_oct_factor* probabilities of measurement success
   for 'square' (four-qubit, intra-cavity) and 'octagon' (eight-qubit, 
   inter-cavity) measurements. Only used when `is_noisy == 1`

Special note: `iterations` will likely have a value which looks low. 
This is because I'm going to batch-parallelize this program by running
~100 jobs each with `iterations` iterations.

1.3 SIMULATION OUTPUT

The simulation result is saved in the out directory in separate files for each lattice size. The directory also contains scripts that visualize the data. For visualization, copy the .txt files to the expt1 directory and type ./get_plot to produce a .ps and .pdf file with the visualizations.
