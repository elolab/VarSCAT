# VaSCAT: The Variants Sequence Context Analysis Toolkit
## Introduction
VarSCAT is an open-source， a command-line based tool written in Python for analyzing variant sequence context. VarSCAT takes a VCF file as input, together with a reference sequence, to give various information about sequence context of normalized variants. **The ambiguous variants analysis module** has functions to give breakpoint ambiguous information about 5’ aligned positions, 3’ aligned positions, effected regions of variants, HGVS nomenclature, distance to 3’ direction adjacent variants and flanking bases of REF and ALT. With a given genomic coordinates, VarSCAT could output sequences of original reference sequence and sequence that contain variants, as well as its complementary sequence. **The tandem repeat region variants analysis module** can analyze sequence context around variants and annotate putative tandem repeat regions that contain variants with user defined parameters for purity, composition, and size of putative tandem repeats.

## Dependencies:
**Currently tested on CentOS 7.9 and Ubuntu 16.04. **<br />
1. Python 3.6.8 (https://www.python.org/), perhaps works on python 2.7.12
2. PyVCF  (https://pyvcf.readthedocs.io/en/latest/,my version is 0.6.8) 
3. Biopython (https://biopython.org/,my version is 1.76, to run the tool, at least version 1.72)
4. Pandas (https://pandas.pydata.org/, my version is 1.1.5)
5. NumPy (https://numpy.org/, my version is 1.19.5)
6. pysam (https://pysam.readthedocs.io/en/latest/index.html, my version 0.18.0) is needed, it may need to be install separately
7. ordered-set (https://github.com/rspeer/ordered-set, 4.0.2)
8. pyfaidx (https://pypi.org/project/pyfaidx/, 0.6.4)


## Usage:
### To get help page of the tool: 
**Main:** python VarSCAT.py -h<br />
**Ambiguous variants analysis module:** python VarSCAT.py -A -h<br />
**Tandam repeat region variants analysis module:** python VarSCAT.py -T -h<br />

## Examples
**Output 5' align positions, 3' align positions, 3' edge positions, HGVS nomenclature, flanking regions of variants, distance to 3' variants**<br />
python VarSCAT.py -A --LRP 1 --HGVS 1 --adjacent 1 --flank 1 --vcf test.vcf.gz --reference chr22.fa --output output<br />
**Output reference sequence and mutated sequence for a specfici location**<br />
python VarSCAT.py -A --mut_seq 1 --location chr22:11318581-11318601 --vcf test.vcf.gz --reference chr22.fa --output output_location<br />
**Output resutls with a bed file**<br />
python VarSCAT.py -A --LRP 1 --HGVS 1 --adjacent 1 --flank 1 --bed regions.bed --vcf test.vcf.gz --reference chr22.fa --output output_bed<br />
**Output 5' align positions, 3' align positions, 3' edge positions and perfect tandem repeat regions** <br />
python VarSCAT.py -A --LRP 1 -T --align_continue 100 --gap_continue 0 --vcf test.vcf.gz --reference chr22.fa --output output_TR
       
