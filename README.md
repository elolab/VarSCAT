# VarConT: The Variants Context Analysis Toolkits
Two modules currently support: complex variants module & repeat regions variants module

To get help page of the tool: python VarCont.py -h
Complex variants module: python VarCont.py -c -h
Repeat regions variants module: python VarCont.py -r -h

Example: python VarCont.py -c --vcf result.vcf --reference reference.fa --location chr19:605124-694960 --based 1 --output_complex my_complex 
         python VarCont.py -r --vcf result.vcf --reference reference.fa --location chr19:605124-694960 --based 1 --output_repeat my_repeat

Currently it only tested and support with Ubuntu 16.04. Window are not supported because it had problem with make. Some dependencies of VarCont can nut run.

Dependencies: 1. PyVCF (my version is 0.6.8)
              2. Biopython (my version is 1.76, to run the tool, at least version 1.72)
              3. Pandas (my version is 0.17.1)
              4. NumPy (my version is 1.16.6)
              
