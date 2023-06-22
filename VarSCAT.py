import sys, getopt
from resource.read_VCF import read_vcf
from resource.read_VCF import complex_module
from resource.read_VCF import repeat_module
import datetime
import pandas as pd

def main():
	vcf_input=''
	reference_f=''
	location=''
	bed=''
	based=1
	output='VarSCAT_output'
	LRP="0"
	HGVS="0"
	flank="0"
	Mut_seq="0"
	Neighbor="0"  
	Complement="0"
	Annotation=''
	rmin=1
	rmax=6 
	mtime=4
	mmper=100.0
	alignment_a=10.0
	M_score=1.0
	MI_score=-1.0
	G_score=-2.0
	R_check=100.0
	S_except=0
	try:
		options, remainder = getopt.gnu_getopt(sys.argv[1:], 'ATh',['Adjacent','TR','help','vcf=','reference=','location=','bed=','based=','output=','LRP=','HGVS=','flank=','mut_seq=','neighbor=','complement=','annotation=','min_unit=','max_unit=','min_time=','min_match_per=','gap_tolerate=','min_score=','match=','mismatch=','gap=','similarity='])		
		if len(sys.argv[1:])==1 and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
			print('''VarSCAT: Variant Sequence Context Annotation Tool (v1.1.0)
Help main:
-A,--Adjacent: adjacent sequence annotation module.
-T,--TR: tandem repeat annotation module.
-h,--help: help page. (-h, -A -h, -T -h)

Two modules can be used together or separate.
If two modules are used together, the commom parameters '--vcf','--reference','--location','--bed','--based' and '--output' should be only announced once. Results of two modules will be merged in one file. If no module is given, the output will be normalized variant list in txt format.

Examples:
# Output 5' align positions, 3' align positions, 3' edge positions, flanking bases of variants, HGVS nomenclature and distance to 3' variants
$python VarSCAT.py -A --LRP 1 --HGVS 1 --flank 1 --neighbor 1 --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output

# Output the reference sequence, the mutated sequence and the reverse complement of mutated sequence for a specfici location
$python VarSCAT.py -A --mut_seq 1 --complement 1 --location chr_test:20-30 --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output_location

# Parse variants for several locations in a bed file
$python VarSCAT.py -A --LRP 1 --HGVS 1 --flank 1 --neighbor 1 --annotation ./data/custom.bed --bed ./data/regions.bed --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output_bed

# Output flanking bases of variants and tandem repeat regions with default setting
$python VarSCAT.py -A --flank 1 -T --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output_TR
''')	
			sys.exit(2)	
		if ('-A' in sys.argv[1:] or '--Adjacent' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
			print('''VarSCAT: Variant Sequence Context Annotation Tool (v1.1.0)
			
Adjacent Sequence Annotation Module:
Required parameters:
--vcf: input VCF file. (The VCF file should be indexed if "--location" or "--bed" is activated, a tbi file of the VCF is required)
--reference: input reference sequencing file. (The reference sequence should be indexed, a fai file is required)
--based: 0-based or 1-based reference coordination. (default:1)
--output: prefix of output file.

Optional parameters:
--location: a genome location needs to be parsed. (format chrx:xxxx-xxxx)
--bed: a bed file contains genome locations need to be parsed.("choromosome", "start", "end" are required)
--LRP: output the 5' aligned (left-most) and 3' aligned (right most) coordinates and 3' edge positions of variants. (default=0, 0:false,1:true)
--HGVS: output the HGVS nomenclature (default=0, 0:false,1:true. Note: According to HGVS recommendation, the reference sequence can only be NCBI Reference Sequence,user should know the corresponding accession and version of the used reference)
--flank: output the flank bases of variants. (default=0, 0:false,1:true)
--neighbor: output the distance to 3' direction nearest variant. (default=0, 0:false,1:true)
--mut_seq: output the reference and mutated sequence based on variants. (default=0, 0:false,1:true. Note: valid with "--location")
--complement: output the reverse complement sequence of mutated sequence. (default=0, 0:false,1:true. Note: valid with "--mut_seq")
--annotation: annotate variants with custom files in bed format. ("choromosome", "start", "end" are required. Additional information can be provided and annotated. Multiple bed files can be used, "custom.bed,custom2.bed". Note: valid with "--location" or "--bed")

-h,--help: help page.
''')
			sys.exit(2)

		if ('-T' in sys.argv[1:] or '--TR' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):	
			print('''VarSCAT: Variant Sequence Context Annotation Tool (v1.1.0)
			
Tandem Repeat Annotation Module:
Required parameters:
--vcf: input VCF file. (The VCF file should be indexed if "--location" or "--bed" is activated, a tbi file of the VCF is required)
--reference: input reference sequencing file. (The reference sequence should be indexed, a fai file is required)
--based: 0-based or 1-based reference coordination. (default:1)
--output: prefix of output file.

Optional parameters:
--location: a genome location needs to be parsed. (format chrx:xxxx-xxxx)
--bed: a bed file contains genome locations need to be parsed.(Three columns: choromosome, start, end)

Advanced parameters:
--min_unit: the minimun size of tandem repeat motifs. (default=1)
--max_unit: the maximum size of tandem repeat motifs. (default=6, larger size will increase the running time)
--min_time: the minimun copy number to call a tandem repeat region. (default=4) 
--match: the match score for motifs aligned with a potential tandem repeat region. (default=1)
--mismatch: the mismatch score for motifs aligned with a potential tandem repeat region. (default=-1)
--gap: the gap penalty for for motifs aligned with a potential tandem repeat region. (default=-2)
--similarity: the minimum similarity between potential repeat units. (default=100, means 100% similarity)
--gap_tolerate: the maximum tolerated gap size (bp) between potential repeat units. (default=0, set -1 for maximum gap of motif size)
--min_score: the minimum alignment sum score for a tandem repeat region. (default=10, set according "--match","--mismatch","--gap")
--min_match_per: the minimum match percentage for a tandem repeat region. (default=100, means 100% of matches)

-h,--help: help page.
''')
			sys.exit(2)
		print("VarSCAT: Variant Sequence Context Annotation Tool (v1.1.0)"+"\n")
		for opt,arg in options:
			if opt in ('--reference'):
				reference_f=arg	
			elif opt in ('--location'):
				location=arg
			elif opt in ('--bed'):
				bed=arg
			elif opt in ('--vcf'):
				vcf_input=arg	
			elif opt in ('--based'):
				based=int(arg)
			elif opt in ('--output'):
				output=arg
		if ('-A' in sys.argv[1:] or '--Adjacent' in sys.argv[1:]):
			for opt,arg in options:
				if opt in ('--LRP'):
					LRP=arg
				elif opt in ('--HGVS'):
					HGVS=arg
				elif opt in ('--flank'):
					flank=arg
				elif opt in ('--mut_seq'):
					Mut_seq=arg
				elif opt in ('--neighbor'):
					Neighbor=arg
				elif opt in ('--complement'):
					Complement=arg
				elif opt in ('--annotation'):
					Annotation=arg	
		if ('-T' in sys.argv[1:] or '--TR' in sys.argv[1:]):
			for opt,arg in options:
				if opt in ('--min_unit'):
					rmin=int(arg)
				elif opt in ('--max_unit'):
					rmax=int(arg)
				elif opt in ('--min_time'):
					mtime=int(arg)
				elif opt in ('--match'):
					M_score=float(arg)
				elif opt in ('--mismatch'):
					MI_score=float(arg)
				elif opt in ('--gap'):
					G_score=float(arg)
				elif opt in ('--similarity'):
					R_check=float(arg)
				elif opt in ('--gap_tolerate'):
					S_except=int(arg)
				elif opt in ('--min_score'):
					alignment_a=float(arg)
				elif opt in ('--min_match_per'):
					mmper=float(arg)
		
		print("Parameters:")
		print("--reference "+reference_f)
		print("--location "+location)
		print("--bed "+bed)
		print("--vcf "+vcf_input)
		print("--based "+str(based))
		print("--output "+output+"\n")
		if ('-A' in sys.argv[1:] or '--Adjacent' in sys.argv[1:]):
			print("--Adjacent")
			print("--LRP "+LRP)
			print("--HGVS "+HGVS)
			print("--flank "+flank)
			print("--neighbor "+Neighbor)
			print("--mut_seq "+Mut_seq)
			print("--complement "+Complement)
			print("--annotation "+Annotation+"\n")
		if ('-T' in sys.argv[1:] or '--TR' in sys.argv[1:]):
			print("--TR")
			print("--min_unit "+str(rmin))
			print("--max_unit "+str(rmax)) 
			print("--min_time "+str(mtime))
			print("--match "+str(M_score))
			print("--mismatch "+str(MI_score))
			print("--gap "+str(G_score))
			print("--similarity "+str(R_check))	
			print("--gap_tolerate "+str(S_except))
			print("--min_score "+str(alignment_a))
			print("--min_match_per "+str(mmper)+"\n")
		if location=="" and Mut_seq=="1":
			print("ERROR: '--mut_seq' needs a valid '--location'")
			sys.exit(2)
		if Complement=="1" and Mut_seq=="0":
			print("ERROR: '--complement' needs a valid '--mut_seq'")
			sys.exit(2)
		if Annotation!="" and (location=="" and bed==""):
			print("ERROR: '--annotation' needs a valid '--location' or '--bed'")
			sys.exit(2)
		if location!="" and bed!="":
			print("ERROR: '--location' and '--bed' can't be valid at the same time")
			sys.exit(2)
		
		start_time = datetime.datetime.now()
		print("Program start at: "+str(start_time)+"\n")
		print("Start Normalization")
		pre_df_set = read_vcf(vcf_input,location,bed,reference_f,based,output)
		print("Normalization Ready"+"\n")
		variant_df = pre_df_set[0]
		sample_list = pre_df_set[1]
		genotype_df = pre_df_set[2]
		lrp_df = pre_df_set[3]
		variant_df_A = []
		variant_df_T = []

		if ('-A' in sys.argv[1:] or '--Adjacent' in sys.argv[1:]):
			print("Start Adjacent Sequence Annotation Module")	
			results_A = complex_module(variant_df,lrp_df,reference_f,based,LRP,HGVS,flank,location,bed,Mut_seq,Neighbor,Complement,Annotation,sample_list,output)
			variant_df_A = results_A[0]
			header_list_A= results_A[1]
			print("Adjacent Sequence Annotation Module Ready"+"\n")	

		if ('-T' in sys.argv[1:] or '--TR' in sys.argv[1:]):
			print("Start Tandem Repeat Annotation Module")
			results_T = repeat_module(variant_df,lrp_df,reference_f,rmin,rmax,mtime,mmper,alignment_a,based,M_score,MI_score,G_score,R_check,S_except)
			variant_df_T = results_T[0]
			header_list_T= results_T[1]
			print("Tandem Repeat Annotation Module Ready"+"\n")	
		
		variant_df=[x + y for x, y in zip(variant_df, genotype_df)]
		if LRP=="1":
			variant_df=list(map(list, zip(*variant_df)))
			variant_df=variant_df+lrp_df
			variant_df=list(map(list, zip(*variant_df)))
			header_list = ["Chromosome","Position","REF","ALT"]+sample_list+["5'_aligned","3'_aligned","3'_edge"]
		else:
			header_list = ["Chromosome","Position","REF","ALT"]+sample_list	

		if len(variant_df_A)!=0:
			variant_df=list(map(list, zip(*variant_df)))
			variant_df=variant_df+variant_df_A
			variant_df=list(map(list, zip(*variant_df)))
			header_list=header_list+header_list_A
		if len(variant_df_T)!=0:
			variant_df=list(map(list, zip(*variant_df)))
			variant_df=variant_df+variant_df_T
			variant_df=list(map(list, zip(*variant_df)))
			header_list=header_list+header_list_T
        
		variant_df = pd.DataFrame(variant_df, columns=[header_list])
		output_file = output+'.txt'
		variant_df.to_csv(output_file, index = False, header=True, sep='\t')
		
		end_time = datetime.datetime.now()
		print("Program end at: "+str(end_time)+"\n")
		print("\n")

	except getopt.GetoptError:
		print('''ERROR: is your input format correct?
		
VarSCAT: Variant Sequence Context Annotation Tool (v1.1.0)
Help main:
-A,--Adjacent: adjacent sequence annotation module.
-T,--TR: tandem repeat annotation module.
-h,--help: help page. (-h, -A -h, -T -h)

Two modules can be used together or separate.
If two modules are used together, the commom parameters '--vcf','--reference','--location','--bed','--based' and '--output' should be only announced once. Results of two modules will be merged in one file. If no module is given, the output will be normalized variant list in txt format.

Examples:
Output 5' align positions, 3' align positions, 3' edge positions, flanking bases of variants, HGVS nomenclature and distance to 3' variants
$python VarSCAT.py -A --LRP 1 --HGVS 1 --flank 1 --neighbor 1 --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output

Output the reference sequence, the mutated sequence and the reverse complement of mutated sequence for a specfici location
$python VarSCAT.py -A --mut_seq 1 --complement 1 --location chr_test:20-30 --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output_location

Parse variants for several locations in a bed file
$python VarSCAT.py -A --LRP 1 --HGVS 1 --flank 1 --neighbor 1 --bed ./data/regions.bed --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output_bed

Output flanking bases of variants and tandem repeat regions with default setting
$python VarSCAT.py -A --flank 1 -T --vcf ./data/test.vcf.gz --reference ./data/test.fa --output output_TR
''')	
		sys.exit(2)

if __name__ == "__main__":
	main()

