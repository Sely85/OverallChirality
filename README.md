OverallChirality
================


This code will compute an overall chirality for your molecule using the Osipov chirality index
(Osipov et al., Mol. Phys., 84, 1193-1206, 1995).
The calculation will take into account all possible permutations of all possible sets of 4 atoms
in your xyz. The output chirality value will be divided for the total number of atoms.

Usage:
        OverallChirality <xyz>


BUILD (Linux)
-------------
Compile the code with OpenMP flag to run it in parallel:

g++ -fopenmp -lm overall_chirality.cpp -o OverallChirality


USAGE
-----
./OverallChirality <xyz>

The code will compute an overall chirality index for your molecule.



EXAMPLE
-------
Execute:
	./OverallChirality limoneneS.xyz
	./OverallChirality limoneneL.xyz

Output:
The resulting chirality index will be -0.0239623 for R-limonene and 0.0239623 for L-limonene.


