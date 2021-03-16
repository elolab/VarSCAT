import sys, getopt
import read_VCF
import read_reference

def reference_fetch(reference_file):
	ref_seq=read_reference.read_reference_sequence(reference_file)
	return ref_seq

def main():
	vcf_input=''
	reference_f=''
	location=''
	based=1
	output_ambiguity='output_ambiguity'
	output_satellite='output_satellite'
	HGVS="0"
	flank="0"
	Mut_seq="0"
	Adjacent="0"  
	Complement="0"
	Transcript="0"
	Translate="0"
	table_def=""
	Functional=""
	rmin=1
	rmax=21 
	esize=5
	alignment_score=3
	gap_score=10
	try:
		options, remainder = getopt.gnu_getopt(sys.argv[1:], 'ash',['ambiguity','satellite','help','vcf=', 'reference=','location=','based=','output_ambiguity=','output_satellite=','HGVS=','flank=','mut_seq=','adjacent= ','complement=','transcript=','translate=','translate_table=','functional=','min_unit=','max_unit=','psize=','alignment_score=','gap_score='])
		if len(sys.argv[1:])==1 and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
			print('''Help main:
VarCont v1.0
-a,--ambiguity: ambigious variants analysis module
-s,--satellite: satellite DNA region variants analysis module
-h,--help: help page

Two modules can be used together, the commom parameters '--vcf','--reference' and '--location' should be only announced once.
''')	
			sys.exit(2)
		
		if ('-a' in sys.argv[1:] or '--ambiguity' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
			print('''Help page of ambigious variants analysis module:
Required parameters:
--vcf: input VCF file 
--reference: input reference sequencing file
--based: 0-based or 1-based reference coordination (default:1)

Optional parameters:
--location: genome location need to be ,parsed (the VCF file should be indexed if -L is activated, tbi file of VCF is required, empty value: annotate all variants in VCF)
--output_ambiguity: ambiguity variants module output prefix
--HGVS: output the HGVS nomenclature (default=0, equal to False. Note: According to HGVS recommendation, the reference sequence can only be NCBI Reference Sequence,user should know the corresponding      
        accession and version of the used reference)
--flank: output the flank bases of variants (default=0, equal to False)
--mut_seq: output the reference and mutated sequence based on variants (default=0,0:off,1:on. Note: valid with --location)
--adjacent: output the flag of adjacent variants (default=0,equal to False)
--complement: output the reverse complement sequence of mutated sequence (default=0, Note: valid with --mut_seq)
--transcript:output the transcription of mutated sequence (default=0, Note: Only valid with --mut_seq, if --complement valid, ouput transcript of reverse complement sequence of mutated sequence)
--translate:output the translation of mutated sequence (default=0, Note: Only valid with --mut_seq,if --complement valid, ouput translarion of reverse complement sequence of mutated sequence)
--translate_table: the translation table (default=1, Standard table. More info:https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11. Note: Only valid --translation is given)
--functional: output the BED file of ambigous variants regions (default=0,0:off,1:on. If --location valid, also output region selected vcf file).

--help: help page
''')
			sys.exit(2)

		if ('-s' in sys.argv[1:] or '--satellite' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):	
			print('''Help page of satellite DNA region variants analysis module:
Required parameters:
--vcf: input VCF file 
--reference: input reference sequencing file
--based: 0-based or 1-based reference coordination (default:1)

Optional parameters:
--location: genome location need to be parsed (the VCF file should be indexed if -L is activated, tbi file of VCF is required, empty value: annotate all variants in VCF)
--output_satellite: satellite DNA variants module output prefix

Advanced parameters:
--min_unit: the minimun size of satellite DNA pattern unit (default=1)
--max_unit: the maximum size of satellite DNA pattern unit (default=20, larger size will increase the running time)
--psize: the maximum size of satellite DNA pattern for perfect match (default=5, the value should between rmin and rmax)
--alignment_score: the minimum average alignment score of each base in satellite DNA regions (default=3, maximum 5, hard-coded minimum=1)
--gap_score: the maximum average gap of each satellite DNA unit (default=1)
--functional: output the BED of satellite DNA variants regions (default=0,0:off,1:on. If --location valid, also output region selected vcf file).

--help: help page
''')
			sys.exit(2)
		
		for opt,arg in options:
			if opt in ('--reference'):
				reference_f=arg	
			elif opt in ('--location'):
				location=arg
			elif opt in ('--vcf'):
				vcf_input=arg	
			elif opt in ('--based'):
				based=int(arg)
			elif opt in ('--functional'):
				Functional=arg
		pre_df = read_VCF.read_vcf(vcf_input,location,Functional,reference_f,based)
		print ("Normalized VCF ready")

		if ('-a' in sys.argv[1:] or '--ambiguity' in sys.argv[1:]):
			print ("processing ambigious variants analysis module")
			for opt,arg in options:
				if opt in ('--output_ambiguity'):
					output_ambiguity=arg
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
				elif opt in ('--transcript'):
					Transcript=arg
				elif opt in ('--translate'):
					Translate=arg	
				elif opt in ('--translate_table'):
					table_def=arg
				
			complex_result = read_VCF.complex_module(pre_df,reference_f,based,HGVS,flank,location,Mut_seq,Adjacent,Complement,Transcript,Translate,table_def,Functional,vcf_input)
			output_file_ambiguity = output_ambiguity+'.csv'
			complex_result.to_csv(output_file_ambiguity, index = False, header=True)
			print ("Ambigious variants analysis module finished")

		if ('-s' in sys.argv[1:] or '--satellite' in sys.argv[1:]):
			print ("processing satellite DNA region variants analysis module")
			for opt,arg in options:
				if opt in ('--output_satellite'):
					output_satellite=arg
				elif opt in ('--min_unit'):
					rmin=int(arg)
				elif opt in ('--max_unit'):
					rmax=int(arg)+1
				elif opt in ('--psize'):
					esize=int(arg)
				elif opt in ('--alignment_score'):
					alignment_score=float(arg)
				elif opt in ('--gap_score'):
					gap_score=float(arg)*10
			
			repeat_result = read_VCF.repeat_module(pre_df,reference_f,rmin,rmax,esize,alignment_score,gap_score,Functional,vcf_input,location,based)
			output_file_satellite = output_satellite+'.csv'
			repeat_result.to_csv (output_file_satellite, index = False, header=True)
			print ("Satellite DNA region variants analysis module finished")

	except getopt.GetoptError:
		print('''Getopt error! help:
VarCont v1.0
-a,--ambiguity: ambigious variants analysis module
-s,--satellite: satellite DNA region variants analysis module
-h,--help: help page
''')
		sys.exit(2)



if __name__ == "__main__":
	main()
















































