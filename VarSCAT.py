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
	Adjacent="0"  
	Complement="0"
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
		options, remainder = getopt.gnu_getopt(sys.argv[1:], 'ATh',['Ambiguity','TR','help','vcf=','reference=','location=','bed=','based=','output=','LRP=','HGVS=','flank=','mut_seq=','adjacent=','complement=','min_unit=','max_unit=','min_time=','min_match_per=','gap_tolerate=','min_score=','match=','mismatch=','gap=','similarity='])		
		if len(sys.argv[1:])==1 and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
			print('''VarSCAT: Variant Sequence Context Annotation Tool (v1.0.0)
Help main:
-A,--Ambiguity: ambigious variant annotation module
-T,--TR: tandem repeat region variant annotation module
-h,--help: help page

Two modules can be used together or separate
If two modules are used together, the commom parameters '--vcf','--reference','--location','--bed','--based' and '--output' should be only announced once. Results of two modules will be merged in one file. If no module is given, the output will be normalized variant list in txt format
''')	
			sys.exit(2)	
		if ('-A' in sys.argv[1:] or '--Ambiguity' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
			print('''VarSCAT: Variant Sequence Context Annotation Tool (v1.0.0)
Ambigious variant annotation module:
Required parameters:
--vcf: input VCF file 
--reference: input reference sequencing file
--based: 0-based or 1-based reference coordination (default:1)
--output: prefix of output file

Optional parameters:
--location: a genome location (format chrx:xxxx-xxxx) need to be parsed. (the VCF file should be indexed if --location is activated, a tbi file of the VCF is required)
--bed: a bed file contains genome locations need to be parsed.(Three columns: choromosome, start, end)
--LRP: output the 5' aligned (left-most) and 3' aligned (right most) coordinates and 3' edge positions of variants. (default=0,equal to False)
--HGVS: output the HGVS nomenclature (default=0, equal to False. Note: According to HGVS recommendation, the reference sequence can only be NCBI Reference Sequence,user should know the corresponding accession and version of the used reference)
--flank: output the flank bases of variants. (default=0, equal to False)
--adjacent: output the distance to 3' direction nearest variant. (Integrated VCF is not supported,default=0,equal to False)
--mut_seq: output the reference and mutated sequence based on variants. (Integrated VCF is not supported,default=0,0:off,1:on. Note: valid with --location)
--complement: output the reverse complement sequence of mutated sequence. (Integrated VCF is not supported,default=0, Note: valid with --mut_seq)

-h,--help: help page
''')
			sys.exit(2)

		if ('-T' in sys.argv[1:] or '--TR' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):	
			print('''VarSCAT: Variant Sequence Context Annotation Tool (v1.0.0)
Tandem repeat region variant annotation module:
Required parameters:
--vcf: input VCF file 
--reference: input reference sequencing file
--based: 0-based or 1-based reference coordination (default:1)
--output: prefix of output file

Optional parameters:
--location: a genome location (format chrx:xxxx-xxxx) need to be parsed. (the VCF file should be indexed if --location is activated, a tbi file of the VCF is required)
--bed: a bed file contains genome locations need to be parsed.(Three columns: choromosome, start, end)

Advanced parameters:
--min_unit: the minimun size of tandem repeat motifs (default=1)
--max_unit: the maximum size of tandem repeat motifs (default=6, larger size will increase the running time)
--min_time: the minimun copy number to call a tandem repeat region (default=4) 
--match: the match score for motifs aligned with a potential tandem repeat region (default=1)
--mismatch: the mismatch score for motifs aligned with a potential tandem repeat region (default=-1)
--gap: the gap penalty for for motifs aligned with a potential tandem repeat region (default=-2)
--similarity: the minimum similarity between potential repeat units (default=100, means 100% similarity)
--gap_tolerate: the maximum tolerated gap size (bp) between potential repeat units (default=0, set -1 for maximum gap of (motif size -1)bp)
--min_score: the minimum alignment sum score for a tandem repeat region (default=10, set according "--match","--mismatch","-gap")
--min_match_per: the minimum match percentage for a tandem repeat region (default=100, means 100% of matches)

-h,--help: help page
''')
			sys.exit(2)
		
		start_time = datetime.datetime.now()
		print("VarSCAT: Variant Sequence Context Annotation Tool"+"\n")
		print("Program start at: "+str(start_time)+"\n")
		print("Parameters:")
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
		print("--reference "+reference_f)
		print("--location "+location)
		print("--bed "+bed)
		print("--vcf "+vcf_input)
		print("--based "+str(based))
		print("--output "+output+"\n")
		print("Start Normalization")
		pre_df_set = read_vcf(vcf_input,location,bed,reference_f,based,output)
		print("Normalization Ready"+"\n")
		variant_df = pre_df_set[0]
		sample_list= pre_df_set[1]
		genotype_df= pre_df_set[2]
		variant_df_A=[]
		variant_df_T=[]

		if ('-A' in sys.argv[1:] or '--Ambiguity' in sys.argv[1:]):
			print("Start Ambiguity Module")
			for opt,arg in options:
				if opt in ('--LRP'):
					LRP=arg
				elif opt in ('--HGVS'):
					HGVS=arg
				elif opt in ('--flank'):
					flank=arg
				elif opt in ('--mut_seq'):
					Mut_seq=arg
				elif opt in ('--adjacent'):
					Adjacent=arg
				elif opt in ('--complement'):
					Complement=arg
			print("--Ambiguity")
			print("--LRP "+LRP)
			print("--HGVS "+HGVS)
			print("--flank "+flank)
			print("--adjacent "+Adjacent)
			print("--mut_seq "+Mut_seq)
			print("--complement "+Complement)	
			results_A = complex_module(variant_df,reference_f,based,LRP,HGVS,flank,location,Mut_seq,Adjacent,Complement,sample_list,output)
			variant_df_A = results_A[0]
			header_list_A= results_A[1]
			print("Ambiguity Module Ready"+"\n")	

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
			print("Start TR module")
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
			print("--min_match_per "+str(mmper))
			results_T = repeat_module(variant_df,reference_f,rmin,rmax,mtime,mmper,alignment_a,based,M_score,MI_score,G_score,R_check,S_except)
			variant_df_T = results_T[0]
			header_list_T= results_T[1]
			print("TR module Ready"+"\n")	

		variant_df=[x + y for x, y in zip(variant_df, genotype_df)]
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
		print('''Getopt error! help:
VarSCAT: Variant Sequence Context Annotation Tool (v1.0.0)
Help main:
-A,--Ambiguity: ambigious variant annotation module
-T,--TR: tandem repeat region variant annotation module
-h,--help: help page

Two modules can be used together or separate.
If two modules are used together, the commom parameters '--vcf','--reference','--location','--bed','--based' and '--output' should be only announced once. Results of two modules will be merged in one file. If no module is given, the output will be normalized variant list in txt format
''')	
		sys.exit(2)

if __name__ == "__main__":
	main()
