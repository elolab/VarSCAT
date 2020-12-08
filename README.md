# VarConT: The Variants Context Analysis Toolkits
Two modules currently support: complex variants module & repeat regions variants module

To get help page of the tool: python VarCont.py -h
Complex variants module: python VarCont.py -c -h
Repeat regions variants module: python VarCont.py -r -h

Currently it only tested and support with Ubuntu 16.04. Window are not supported because it had problem with make. Some dependencies of VarCont can nut run.

Dependencies: 1. PyVCF (my version is 0.6.8)
              2. Biopython (my version is 1.76, to run the tool, at least version 1.72)
              3. Pandas (my version is 0.17.1)
              4. NumPy (my version is 1.16.6)
              