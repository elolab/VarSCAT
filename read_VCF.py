import sys
import pandas as pd
import vcf
import LRP_system
import var_context 
import mut_ref_var
import Complex_variants
import Repeat_context
import HGVS_convert
import read_reference

def Redundant_variant(P,R,A,ref,bas):
	while R[-1]==A[-1]:
		R=R[:len(R)-1]
		A=A[:len(A)-1]
		if len(R)==0:
			P=int(P)-1
			R=ref[P-bas].upper()
			A=ref[P-bas].upper()+A
		elif len(A)==0:
			P=int(P)-1
			R=ref[P-bas].upper()+R
			A=ref[P-bas].upper()

	if len(R)>1 and len(A)>1:
		while R[0]==A[0] and len(R)>1 and len(A)>1:
			R=R[1:]
			A=A[1:]
			P=str(int(P)+1)
	return P,R,A

def create_variant_list(vcf_reader_f):
	variant_list_f = []
	for i in vcf_reader_f:
		chromosome_f=str(i.CHROM)
		position_f = str(i.POS)
		REF_f = str(i.REF).upper()
		ALT_f = str(i.ALT).upper()
		ALT_f = ALT_f.replace("[","")
		ALT_f = ALT_f.replace("]","")
		ALT_f = ALT_f.replace(" ","")
		
		sample_name_list_f = []
		genotype_list_f = []
		if len(vcf_reader_f.samples)!=0:
			for s in vcf_reader_f.samples:
				sample_name_f = s
				sample_name_list_f.append(sample_name_f)
				genotype_s_f = str(i.genotype(sample_name_f)['GT'])
				genotype_list_f.append(genotype_s_f)
		else:
			sample_name_f = "Genotype"
			sample_name_list_f.append(sample_name_f)
			genotype_s_f = "0/0"
			genotype_list_f.append(genotype_s_f)
		if "," in ALT_f:
			ALT_multi = ALT_f.split(",")
			for j in ALT_multi:
				ALT_2 = j
				variant_info1 = [chromosome_f, position_f, REF_f, ALT_2]
				variant_info = variant_info1 + genotype_list_f
				variant_list_f.append(variant_info)
		else:
			variant_info1 = [chromosome_f, position_f, REF_f, ALT_f]
			variant_info = variant_info1 + genotype_list_f
			variant_list_f.append(variant_info)
	return variant_list_f,sample_name_list_f

def read_vcf(vcf_file,location,bed,ref_file,based,output):
	if len(location) == 0 and len(bed) == 0:
		try:
			vcf_reader = vcf.Reader(filename=vcf_file)
		except ValueError: 
			print ("Error1:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
			sys.exit(2)
		parse_vcf=create_variant_list(vcf_reader)
		variant_list=parse_vcf[0]
		sample_name_list=parse_vcf[1]
		
	elif len(location) != 0 and len(bed) == 0:	
		chromosome_read = location.split(':')[0]
		start_read = location.split(':')[1].split('-')[0]
		end_read = location.split(':')[1].split('-')[1]
		try:
			vcf_reader = vcf.Reader(filename=vcf_file).fetch(str(chromosome_read), int(start_read), int(end_read))
		except ValueError: 
			print ("Error1:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
			sys.exit(2)
		parse_vcf=create_variant_list(vcf_reader)
		variant_list=parse_vcf[0]
		sample_name_list=parse_vcf[1]

	elif len(location) == 0 and len(bed) != 0:
		bed_file=open(bed,"r")
		bed_file_lines=bed_file.readlines()
		variant_list=[]
		sample_name_list=[]
		for item in bed_file_lines:
			item=item.replace("\n","")
			item=item.split("\t")
			chromosome_read = item[0]
			start_read = item[1]
			end_read = item[2]
			try:
				vcf_reader = vcf.Reader(filename=vcf_file).fetch(str(chromosome_read), int(start_read), int(end_read))
			except ValueError: 
				print ("Error1:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
				sys.exit(2)
			parse_vcf=create_variant_list(vcf_reader)
			variant_list=variant_list+parse_vcf[0]
		sample_name_list=parse_vcf[1]
		bed_file.close()
			
	variant_list_df = []
	chr_name=""
	for i in range(0,len(variant_list)):
		chromosome=variant_list[i][0]
		if chr_name!=chromosome:
			chr_name=chromosome
			reference=read_reference.read_reference_sequence(ref_file,chr_name)
		if len(variant_list[i][2])!=1 or len(variant_list[i][3])!=1:
			Normalized_variant = Redundant_variant(variant_list[i][1],variant_list[i][2],variant_list[i][3],reference,based)
			A=Normalized_variant[0]
			B=Normalized_variant[1]
			C=Normalized_variant[2]
			item = {"Chromosome":variant_list[i][0],"Position":A,"REF":B,"ALT":C}
			for j in range(4,len(variant_list[i])):
				item[sample_name_list[j-4]]=variant_list[i][j]
			variant_list_df.append(item)

		else:
			item = {"Chromosome":variant_list[i][0],"Position":variant_list[i][1],"REF":variant_list[i][2],"ALT":variant_list[i][3]}
			for j in range(4,len(variant_list[i])):
				item[sample_name_list[j-4]]=variant_list[i][j]
			variant_list_df.append(item)

	variant_list_df = pd.DataFrame(variant_list_df)
	variant_list_df = variant_list_df[["Chromosome","Position","REF","ALT"]+sample_name_list]
	
	return variant_list_df,sample_name_list

def complex_module(df,ref,bas,lrp,HGVS,fla,loc,Mut,adj,com,vcf,sample,output):
	item_LRP = LRP_system.LRP(df,ref,bas)
	df["5'_aligned"] = item_LRP[0]
	df["3'_aligned"] = item_LRP[1]
	df["3'_edge"] = item_LRP[2]

	if fla=="1":
		item_varcon = var_context.mutation_seq(df,ref,bas)
		df["ref_sequence"] = item_varcon[0]
		df["mut_sequence"] = item_varcon[1]
	
	if loc!="" and Mut=="1" and len(sample)==1:
		item_mutref = mut_ref_var.mut_ref(df,loc,ref,com,bas,output)

	if HGVS=="1":
		item_HGVS=HGVS_convert.HGVS_transform(df,ref,bas)
		df["HGVS"] = item_HGVS
	
	if adj=="1":
		item_complex = Complex_variants.combine_complex(df)
		df["distance_3_nearest_Var(bp)"] = item_complex
	
	if lrp=="0":
		del df["5'_aligned"]
		del df["3'_edge"]
		del df["3'_aligned"]

	if loc=="" and Mut=="1":
		print("Non-valid --location is not supported for --mut_seq")
	if len(sample)>1 and Mut=="1":
		print("Integrated VCF is not supported for --mut_seq")
			
	return df

def repeat_module(df,ref,min_size,max_size,min_time,min_per,alignAS,vcf,loc,bas,Mscore,MIscore,OGscore,Rcheck,output,Sexcept):
	min_per=min_per/100
	Rcheck=Rcheck/100
	item_repeat = Repeat_context.search_repeats(df,ref,min_size,max_size,min_time,min_per,alignAS,bas,Mscore,MIscore,OGscore,Rcheck,Sexcept)
	df["Motifs"] =item_repeat[0]
	df["Copy_number"] =item_repeat[1]
	df["Size"] =item_repeat[2]
	df["Start"] =item_repeat[3]
	df["End"] =item_repeat[4]
	df["Repeat_Score"] =item_repeat[5]
	df["Alignment_Score"] =item_repeat[6]
	df["Match%"] =item_repeat[7]
	df["Mismatch%"] =item_repeat[8]
	df["Gap%"] =item_repeat[9]
	df["Repeat_GC%"] =item_repeat[10]

	return df

