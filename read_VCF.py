import sys
import pandas as pd
import vcf
import UPS_system
import var_context 
import mut_ref_var
import Complex_variants
import Repeat_context
import function_output
import HGVS_convert
import read_reference

def Redundant_variant(P,R,A,ref,bas):
	while R[-1]==A[-1]:
		R=R[:len(R)-1]
		A=A[:len(A)-1]
		if len(R)==0:
			P=int(P)-1
			R=ref[P-bas]
			A=ref[P-bas]+A
		elif len(A)==0:
			P=int(P)-1
			R=ref[P-bas]+R
			A=ref[P-bas]

	if len(R)>1 and len(A)>1:
		while R[0]==A[0] and len(R)>1 and len(A)>1:
			R=R[1:]
			A=A[1:]
			P=str(int(P)+1)
	return P,R,A

def read_vcf(vcf_file,location,func,ref_file,based):
	if len(location) == 0:
		vcf_reader = vcf.Reader(filename=vcf_file)
	else:	
		chromosome_read = location.split(':')[0]
		start_read = location.split(':')[1].split('-')[0]
		end_read = location.split(':')[1].split('-')[1]
		try:
			vcf_reader = vcf.Reader(filename=vcf_file).fetch(str(chromosome_read), int(start_read), int(end_read))
			if func!="":
				vcf_readerF = vcf.Reader(filename=vcf_file).fetch(str(chromosome_read), int(start_read), int(end_read))
				vcf_writer = vcf.Writer(open('VCF_in_region.vcf', 'w'), vcf_readerF)
				for record in vcf_readerF:
					vcf_writer.write_record(record)
				vcf_writer.close()
		except ValueError: 
			print ("Error1:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")				
	
	sample_name = vcf_reader.samples[0]

	variant_list = []
	for i in vcf_reader:
		chromosome=str(i.CHROM)
		position = str(i.POS)
		REF = str(i.REF)
		ALT = str(i.ALT)
		ALT = ALT.replace("[","")
		ALT = ALT.replace("]","")
		ALT = ALT.replace(" ","")
		genotype = str(i.genotype(sample_name)['GT'])

		if "," in ALT:
			ALT_multi = ALT.split(",")
			for j in ALT_multi:
				ALT_2 = j
				variant_info = [chromosome,position, REF, ALT_2, genotype]
				variant_list.append(variant_info)
		else:
			variant_info = [chromosome,position, REF, ALT, genotype]
			variant_list.append(variant_info)
		
	variant_list_df = []
	chr_name=""
	for i in range(0,len(variant_list)):
		chromosome=variant_list[i][0]
		if chr_name!=chromosome:
			chr_name=chromosome
			reference=read_reference.read_reference_sequence(ref_file,chr_name)

		if len(variant_list[i][2])!=1 and len(variant_list[i][3])!=1:
			Normalized_variant = Redundant_variant(variant_list[i][1],variant_list[i][2],variant_list[i][3],reference,based)
			A=Normalized_variant[0]
			B=Normalized_variant[1]
			C=Normalized_variant[2]
			item = {"Chromosome":variant_list[i][0],"Position":A,"REF":B,"ALT":C,"Genotype":variant_list[i][4]}
			variant_list_df.append(item)

		else:
			item = {"Chromosome":variant_list[i][0],"Position":variant_list[i][1],"REF":variant_list[i][2],"ALT":variant_list[i][3],"Genotype":variant_list[i][4]}
			variant_list_df.append(item)

	variant_list_df = pd.DataFrame(variant_list_df)
	variant_list_df = variant_list_df[["Chromosome","Position","REF","ALT","Genotype"]]

	print("Normalization ready")

	return variant_list_df

def complex_module(df,ref,bas,HGVS,fla,loc,Mut,adj,com,tap,tal,table_def,func,vcf):
	item_UPS = UPS_system.UPS(df,ref,bas)
	df["UPS_L"] = item_UPS[0]
	df["UPS_R"] = item_UPS[1]
	df["Right align position"] = item_UPS[2]
	df["UPS_LEFT"] = item_UPS[3]
	df["UPS_RIGHT"] = item_UPS[4]

	if fla=="1":
		item_varcon = var_context.mutation_seq(df,ref,bas)
		df["ref_sequence"] = item_varcon[0]
		df["mut_sequence"] = item_varcon[1]
	
	if loc!="" and Mut=="1":
		item_mutref = mut_ref_var.mut_ref(df,loc,ref,com,tap,tal,table_def,based)

	if HGVS=="1":
		item_HGVS=HGVS_convert.HGVS_transform(df,ref,bas)
		df["HGVS"] = item_HGVS
	
	if adj=="1":
		item_complex = Complex_variants.combine_complex(df)
		df["Adjacent variants"] = item_complex[0]
		df["distance_3'_nearest_Var(bp)"] = item_complex[1]

	if func!="":
		function_output.Cfunc(df)
	
	del df["UPS_LEFT"]
	del df["UPS_RIGHT"]
			
	return df

def repeat_module(df,ref,min_size,max_size,exact_size,alignAS,gapAS,func,vcf,loc,bas):
	item_repeat = Repeat_context.search_repeats(df,ref,min_size,max_size,exact_size,alignAS,gapAS,bas)

	if func!="":
		function_output.Rfunc(df)
	
	return item_repeat


 






























