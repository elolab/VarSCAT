# VarSCAT: The Variants Sequence Context Analysis Toolkit
## Introduction
VarSCAT is an open-source, a command-line based tool written in Python for analyzing variant sequence context. VarSCAT takes a VCF file as input, together with a reference sequence, to give various information about sequence context of normalized variants. **The ambiguous variants analysis module** has functions to give breakpoint ambiguous information about 5’ aligned positions, 3’ aligned positions, effected regions of variants, HGVS nomenclature, distance to 3’ direction adjacent variants and flanking bases of REF and ALT. With a given genomic coordinates, VarSCAT could output sequences of original reference sequence and sequence that contain variants, as well as its complementary sequence. **The tandem repeat region variants analysis module** can analyze sequence context around variants and annotate putative tandem repeat regions that contain variants with user defined parameters for purity, composition, and size of putative tandem repeats.

## Dependencies:
**Currently tested on CentOS 7.9**<br />
**Install dependencies:** `pip install -r requirements.txt`
1. Python 3.6.8 (https://www.python.org/)
2. PyVCF>=0.6.8 (https://pyvcf.readthedocs.io/en/latest/) 
3. Biopython>=1.72 (https://biopython.org/)
4. Pandas>=1.1.5 (https://pandas.pydata.org/)
5. pysam>=0.18.0 (https://pysam.readthedocs.io/en/latest/index.html)
6. ordered-set>=4.0.2 (https://github.com/rspeer/ordered-set)
7. pyfaidx>=0.6.4 (https://pypi.org/project/pyfaidx/)


## Usage:
### To get help page of the tool: 
**Main:** `python VarSCAT.py -h`<br />
**Ambiguous variants analysis module:** `python VarSCAT.py -A -h`<br />
```
Ambigious variants analysis module:
Required parameters:`
--vcf: input VCF file
--reference: input reference sequencing file
--based: 0-based or 1-based reference coordination (default:1)
--output: prefix of output file
Optional parameters:
--location: a genome location (format chrx:xxxx-xxxx) need to be parsed. (the VCF file should be indexed)
--bed: a bed file contains genome locations need to be parsed.(Three columns: choromosome, start, end)
--LRP: output the 5' aligned (left-most) and 3' aligned (right most) coordinates and 3' edge positions of variants. (default=0,equal to False)
--HGVS: output the HGVS nomenclature (default=0, equal to False)
--flank: output the flank bases of variants. (default=0, equal to False)
--adjacent: output the distance to 3' direction nearest variant. (Integrated VCF is not supported,default=0,equal to False)
--mut_seq: output the reference and mutated sequence based on variants.(Integrated VCF not supporte,default=0,0:off,1:on.Note: valid with --location)
--complement: output the reverse complement sequence of mutated sequence. (Integrated VCF not supporte,default=0, Note: valid with --mut_seq)
-h,--help: help page
```
**Tandam repeat region variants analysis module:** `python VarSCAT.py -T -h`<br />
```
Tandem repeat region variants analysis module:
Required parameters:
--vcf: input VCF file 
--reference: input reference sequencing file
--based: 0-based or 1-based reference coordination (default:1)
--output: name of output file
Optional parameters:
--location: a genome location (format chrx:xxxx-xxxx) need to be parsed. (the VCF file should be indexed)
--bed: a bed file contains genome locations need to be parsed.(Three columns: choromosome, start, end)
Advanced parameters:
--min_unit: the minimun size of tandem repeat pattern unit. (default=1)
--max_unit: the maximum size of tandem repeat pattern unit. (default=6, larger size will increase the running time)
--min_time: the minimun repeat time to call a tandem repeat. (default=4) 
--match: match score for local alignment of potential repeat unit. (default=1)
--mismatch: mismatch score for local alignment of potential repeat unit. (default=-1)
--gap: gap penalty for local alignment of potential repeat unit. (default=-2)
--align_continue: the minimum similarity percentage of a potential repeat unit that allows local alignment algorithm continue. (default=100, means 100% of similarity)
--gap_continue: the maximum tolerated gap size (bp) between potential repeat units that allows local alignment algorithm continue. (default=0, set to -1 for size of current repeat unit -1))
--min_score: the minimum alignment score of a tandem repeat region. (default=10, set according "--match","--mismatch","-gap")
--min_match_per: the minimum match percentage of a tandem repeat region. (default=100, means 100% of matches)
-h,--help: help page
```
### Examples with test files in data folder
**Output 5' align positions, 3' align positions, 3' edge positions, HGVS nomenclature, flanking bases of variants, distance to 3' variants**<br />
`python VarSCAT.py -A --LRP 1 --HGVS 1 --adjacent 1 --flank 1 --vcf ./data/test.vcf.gz --reference ./data/chr22.fa --output output`<br />
**Output the reference sequence, the mutated sequence and the reverse complement of mutated sequence for a specfici location**<br />
`python VarSCAT.py -A --mut_seq 1 --complement 1 --location chr22:11318581-11318601 --vcf ./data/test.vcf.gz --reference ./data/chr22.fa --output output_location`<br />
**Parse variants for several locations in a bed file**<br />
`python VarSCAT.py -A --LRP 1 --HGVS 1 --adjacent 1 --flank 1 --bed ./data/regions.bed --vcf ./data/test.vcf.gz --reference ./data/chr22.fa --output output_bed`<br />
**Output flanking bases of variants and tandem repeat regions with default setting** <br />
`python VarSCAT.py -A --flank 1 -T --vcf ./data/test.vcf.gz --reference ./data/chr22.fa --output output_TR`
       
