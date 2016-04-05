# Quantum-Toolbox
This are some routines I wrote for my bachelor thesis. They are not intended to be used by anyone else, and are not fully organized and optimized yet.

---Description---

quit_toolbox:
basic routines:
matrix multiplication , logarithm, square roots, trace, hermitian conjugates. 
generation of random states, 
relative entropies , Hilbert schmidt norm, negativity  ect. ect. 

parametrized_state:
returns a fully parametrized mixed state separable relative to a split, or a product state or a simple pure state.

quit_unitary:
implementation of the Jarlskog parametrization of the unitary group U(n)

opt-funktion:
implementation of quantum discord, relative entropy of entaglement, information deficit, geometric discord for entanglement distribution with qubit sized ancilla.

quit_import:
reads a txtfile with a mathematica exported matrix and converts it into an array.

---this code was completely written by myself ---

compile with gcc.4.9

compiler flags: -lstdc++ -lm -lopenblas -lgsl -llapacke -std=c++11
