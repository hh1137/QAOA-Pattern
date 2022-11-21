# QAOA-Pattern generator and correctness check
This generator produces the pattern of compiling bipartite all-to-all interactions at 2xUnits Google Sycamore architectures.

## Dependencies
Building this software requires `Qiskit`, Version `0.34.2`

## Usage and parameters
To use the QAOA-Pattern generator, use the `python QAOAcheck.py` or `python QAOAcheck.py --number 4 --output True --check True`

All the parameter are followings:

'--number n' indicates the input QAOA problem graph is a bipartite clique graph with 2n vertices. For example you can use 'python QAOAcheck.py --number 5` to inidicate a QAOA-10 bipartite clique program graph. n has to be larger than 2. 


'--check' invokes the correctness check of the program if this set as `True`, the default is `False`.

'--output' invokes the output qasm file if this set as `True`, the default is `False`.

