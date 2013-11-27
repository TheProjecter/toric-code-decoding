----------------------------------
      TORIC & [[4,2,2]]-CONCATENATED TORIC CODE SIMULATOR

	Tobias Hoelzer
	September 2013
	Institute for Quantum Information
	RWTH Aachen University
	Bachelor Thesis: Study of [[4,2,2]]-concatenated toric code for a scalable circuit-QED architecture


	based on the works of

      Martin Suchara
      UC Berkeley, 2012
----------------------------------

WHAT IS IT DOING?

All programs determine the failure rates versus error probability for different lattice sizes for a specific amount of iterations and save them to files.
Every program features a specific architecture.


VERSION DESCRIPTION:

The threshold simulations vary in the architecture and the syndromes:

The code is either toric or [[4,2,2]]-concatenated.
The syndrome measurements are either perfect or noisy.



FEATURES:

The toric codes do not differ a lot from Martin Suchara's Code but only consider X-errors with a gauge qubit.
All codes only consider X-errors with a gauge qubit.
All codes feature multithreading.



HOW TO GET STARTED:

The toric codes are easier to understand, but are not very well commented.
The concatenated codes are harder programmingwise, but feature detailed comments.

It may help to read chapter 3 Numerics of the bachelor thesis.



AND LAST, BUT NOT LEAST:

May the force be with you!


