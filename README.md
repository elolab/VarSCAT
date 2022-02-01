# VaSCAT: The Variants Sequence Context Analysis Toolkit
Two modules currently support: complex variants module & repeat regions variants module

To get help page of the tool: 
1. python VarCont.py -h
2. Ambiguous variants module: python VarCont.py -A -h
3. Tandam repeat variants module: python VarCont.py -T -h


Example: python VarSCAT.py -A --vcf result.vcf --reference reference.fa --location chr19:605124-694960 --output_complex my_complex 

Currently it only tested and support with Ubuntu 16.04. Window are not supported because it had problem with make. Some dependencies of VarSCAT can nut run.

Dependencies: 
1. PyVCF (my version is 0.6.8)
2. Biopython (my version is 1.76, to run the tool, at least version 1.72)
3. Pandas (my version is 0.17.1)
4. NumPy (my version is 1.16.6)
5. pysam is needed, it may need to be install separately              
