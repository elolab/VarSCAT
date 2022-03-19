# VaSCAT: The Variants Sequence Context Analysis Toolkit
Two modules currently support: Ambiguous variants analysis module & repeat regions variants analysis module

To get help page of the tool: 
1. python VarSCAT.py -h
2. Ambiguous variants analysis module: python VarSCAT.py -A -h
3. Tandam repeat variants analysis module: python VarSCAT.py -T -h
 
Currently it only tested and support with Ubuntu 16.04. Window are not supported because it had problem with make. Some dependencies of VarSCAT can nut run.

Currently tested on CentOS 7.9 and Ubuntu 16.04

Dependencies:
1. Python 3.6.8, perhaps works on python 2.7.12
2. PyVCF (my version is 0.6.8) 
3. Biopython (my version is 1.76, to run the tool, at least version 1.72)
4. Pandas (my version is 1.1.5)
5. NumPy (my version is 1.19.5)
6. pysam (my version 0.18.0) is needed, it may need to be install separately
7. ordered-set (4.0.2)
8. pyfaidx (0.6.4)
       
