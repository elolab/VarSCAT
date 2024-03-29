import sys
import vcf
from .LRP_function import LRP
from .var_context import mutation_seq
from .mut_ref_var import mut_ref
from .Complex_variants import combine_complex
from .Repeat_context import search_repeats
from .HGVS_convert import HGVS_transform
from .read_reference import read_reference_sequence
from .anno_custom import custom_anno

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
	genotype_list_f = []
	
	sample_name_list_f = []
	if len(vcf_reader_f.samples)!=0:
		for s in vcf_reader_f.samples:
			sample_name_list_f.append(s)
	else:
		sample_name_list_f.append("Genotype")
	sample_name_list_f.insert(0,"ID")
	
	for i in vcf_reader_f:
		chromosome_f=str(i.CHROM)
		position_f = str(i.POS)
		REF_f = str(i.REF).upper()
		ALT_f = str(i.ALT).upper()
		ALT_f = ALT_f.replace("[","")
		ALT_f = ALT_f.replace("]","")
		ALT_f = ALT_f.replace(" ","")
		ID_f = str(i.ID)
		
		genotype_list_f2 = []
		if len(vcf_reader_f.samples)!=0:
			for s in vcf_reader_f.samples:
				genotype_s_f = str(i.genotype(s)['GT'])
				genotype_list_f2.append(genotype_s_f)
		else:
			genotype_s_f = "0/0"
			genotype_list_f2.append(genotype_s_f)
			
		if "," in ALT_f:
			ALT_multi = ALT_f.split(",")
			ID_multi = ID_f.split(",")
			for j in range(0,len(ALT_multi)):
				ALT_2 = ALT_multi[j]
				if len(ID_multi)==len(ALT_multi):
					ID_2 = ID_multi[j]
				else:
					ID_2 = "."
				variant_info1 = [chromosome_f, position_f, REF_f, ALT_2]
				variant_list_f.append(variant_info1)
				genotype_list_f3 = genotype_list_f2.copy()
				genotype_list_f3.insert(0,ID_2)
				genotype_list_f.append(genotype_list_f3)
		else:
			variant_info1 = [chromosome_f, position_f, REF_f, ALT_f]
			variant_list_f.append(variant_info1)
			if ID_f == "None":
			    ID_f = "."
			genotype_list_f2.insert(0,ID_f)
			genotype_list_f.append(genotype_list_f2)
			
	return variant_list_f,sample_name_list_f,genotype_list_f

def read_vcf(vcf_file,location,bed,ref_file,based,output):
	if len(location) == 0 and len(bed) == 0:
		try:
			vcf_reader = vcf.Reader(filename=vcf_file)
		except ValueError: 
			print ("Error:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
			sys.exit(2)
		parse_vcf=create_variant_list(vcf_reader)
		variant_list=parse_vcf[0]
		sample_name_list=parse_vcf[1]
		geno_type_list=parse_vcf[2]
		
	elif len(location) != 0 and len(bed) == 0:	
		chromosome_read = location.split(':')[0]
		start_read = location.split(':')[1].split('-')[0]
		end_read = location.split(':')[1].split('-')[1]
		try:
			vcf_reader = vcf.Reader(filename=vcf_file).fetch(str(chromosome_read), int(start_read)-1, int(end_read))
		except ValueError: 
			print ("Error:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
			sys.exit(2)
		parse_vcf=create_variant_list(vcf_reader)
		variant_list=parse_vcf[0]
		sample_name_list=parse_vcf[1]
		geno_type_list=parse_vcf[2]
		
		if len(variant_list)==0:
			print("WARNING: No variant in the selected region "+chromosome_read+":"+start_read+"-"+end_read)

	elif len(location) == 0 and len(bed) != 0:
		bed_file=open(bed,"r")
		bed_file_lines=bed_file.readlines()
		variant_list=[]
		sample_name_list=[]
		geno_type_list=[]
		for item in bed_file_lines:
			item=item.replace("\n","")
			item=item.split("\t")
			chromosome_read = item[0]
			start_read = item[1]
			end_read = item[2]
			try:
				vcf_reader = vcf.Reader(filename=vcf_file).fetch(str(chromosome_read), int(start_read)-1, int(end_read))
			except ValueError: 
				print ("Error:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
				sys.exit(2)
			parse_vcf=create_variant_list(vcf_reader)
			variant_list=variant_list+parse_vcf[0]
			geno_type_list=geno_type_list+parse_vcf[2]
			if len(parse_vcf[0])==0:
				print("WARNING: No variant in the selected region "+chromosome_read+":"+start_read+"-"+end_read)
		sample_name_list=parse_vcf[1]
		bed_file.close()
			
	variant_list_df = []
	chr_name=""
	for i in range(0,len(variant_list)):
		chromosome=variant_list[i][0]
		if chr_name!=chromosome:
			chr_name=chromosome
			reference=read_reference_sequence(ref_file,chr_name)
		if len(variant_list[i][2])!=1 or len(variant_list[i][3])!=1:
			Normalized_variant = Redundant_variant(variant_list[i][1],variant_list[i][2],variant_list[i][3],reference,based)
			A=Normalized_variant[0]
			B=Normalized_variant[1]
			C=Normalized_variant[2]
			item = [variant_list[i][0],A,B,C]
			variant_list_df.append(item)

		else:
			variant_list_df.append(variant_list[i])
	
	item_LRP = LRP(variant_list_df,ref_file,based)
	variant_list_df_lrp=[]
	variant_list_df_lrp.append(item_LRP[0])
	variant_list_df_lrp.append(item_LRP[1])
	variant_list_df_lrp.append(item_LRP[2])

	return variant_list_df,sample_name_list,geno_type_list,variant_list_df_lrp

def complex_module(df,df_lrp,ref,bas,lrp,HGVS,fla,loc,ibed,Mut,nei,com,anno,samples,output):
	if len(samples)>2 and Mut=="1":
		print("WARNING: Input VCF is a sample-integrated VCF, may not suitable for '--mut_seq'")
	
	df2=[]	
	df2_header=[]
	if fla=="1" or Mut=="1" or HGVS=="1" or nei=="1" or anno!="":
		df=list(map(list, zip(*df)))
		df=df+df_lrp
		df=list(map(list, zip(*df)))
	
		if fla=="1":
			item_varcon = mutation_seq(df,ref,bas)
			df2.append(item_varcon[0])
			df2.append(item_varcon[1])
			df2_header.append("ref_sequence")
			df2_header.append("alt_sequence")
		
		if Mut=="1":
			item_mutref = mut_ref(df,loc,ref,com,bas,output)

		if HGVS=="1":
			item_HGVS=HGVS_transform(df,ref,bas)
			df2.append(item_HGVS)
			df2_header.append("HGVS")
		
		if nei=="1":
			item_complex = combine_complex(df)
			df2.append(item_complex)
			df2_header.append("distance_3_nearest_Var(bp)")
		
		if anno!="":
			anno = anno.split(",")
			for i in range(0,len(anno)):
				item_anno = custom_anno(df,loc,ibed,anno[i])
				df2.append(item_anno[0])
				df2.append(item_anno[1])
				if i==0:
					df2_header.append("Ann_loc")
					df2_header.append("Ann_info")
				else:
					anno_header_1 = "Ann_loc_" + str(i)
					df2_header.append(anno_header_1)
					anno_header_2 = "Ann_info_" + str(i)
					df2_header.append(anno_header_2)
	return df2,df2_header

def repeat_module(df,df_lrp,ref,min_size,max_size,min_time,min_per,alignAS,bas,Mscore,MIscore,OGscore,Rcheck,Sexcept):
	min_per=min_per/100
	Rcheck=Rcheck/100
	
	df=list(map(list, zip(*df)))
	df=df+df_lrp
	df=list(map(list, zip(*df)))
	
	item_repeat = search_repeats(df,ref,min_size,max_size,min_time,min_per,alignAS,bas,Mscore,MIscore,OGscore,Rcheck,Sexcept)

	df2=[]
	df2_header=[]
	df2.append(item_repeat[0])
	df2.append(item_repeat[1])
	df2.append(item_repeat[2])
	df2.append(item_repeat[3])
	df2.append(item_repeat[4])
	df2.append(item_repeat[5])
	df2.append(item_repeat[6])
	df2.append(item_repeat[7])
	df2.append(item_repeat[8])
	df2.append(item_repeat[9])
	df2.append(item_repeat[10])
	df2.append(item_repeat[11])

	df2_header.append("Motifs")
	df2_header.append("Copy_number")
	df2_header.append("Size")
	df2_header.append("Start")
	df2_header.append("End")
	df2_header.append("Repeat_Score")
	df2_header.append("Alignment_Score")
	df2_header.append("Match%")
	df2_header.append("Mismatch%")
	df2_header.append("Gap%")
	df2_header.append("Repeat_GC%")
	df2_header.append("Copy_number_change")

	return df2, df2_header

