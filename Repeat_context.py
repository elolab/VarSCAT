import sys
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import Align
from Bio import motifs
import repeatAlign 
import remove_redundancy
import Local_GC
import read_reference

def catching_variant(position_l1, REF_l1, ALT_l1, Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1,based1):
	if REF_l1[0]==ALT_l1[0]:
		if len(REF_l1) > len(ALT_l1):
			variant_start = int(position_l1)+1-based1
			variant_pattern = REF_l1[1:]
			variant_length = len(variant_pattern)
		elif len(REF_l1) < len(ALT_l1):
			variant_start = int(position_l1)-based1
			variant_pattern = Ref1[variant_start:variant_start+1]
			variant_length = len(variant_pattern)
	elif REF_l1[0]!=ALT_l1[0]:
		variant_start = int(position_l1)-based1
		variant_pattern = REF_l1
		variant_length = len(variant_pattern)

	Repeat_raw_L1 = sizing_window(variant_start,variant_pattern,variant_length,Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1,based1) 
	return Repeat_raw_L1

def catching_variant_SNV(position_l1, REF_l1, ALT_l1, Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1,based1):
	variant_start = int(position_l1)-based1
	variant_pattern = REF_l1
	variant_length = len(REF_l1)
	Repeat_raw_L1 = sizing_window(variant_start,variant_pattern,variant_length,Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1,based1) 
	return Repeat_raw_L1

def sizing_window(start,pattern,length, Ref2,min_size2,max_size2,exact_size2,alignAS2,gapAS2,based2):
	Repeat_df_L2 = pd.DataFrame(columns=["Motifs", "Period","RepeatS","Size","Start","End","Position_O","AlignmentS","GapS"])
	for i in range(min_size2,max_size2):  
		if i <=exact_size2: 
			Repeat_raw_L2 = sliding_window_micro(i,start,length,Ref2,based2)
			if Repeat_raw_L2.shape[0]>1:
				Repeat_raw_L2_sort = remove_redundancy.clean_df(Repeat_raw_L2)
				Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2_sort,ignore_index = True)
			elif Repeat_raw_L2.shape[0]==1:
				Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2,ignore_index = True)
		elif i >exact_size2:
			Repeat_raw_L2 = sliding_window_mini(i,start,length,Ref2,alignAS2,gapAS2,based2)
			if Repeat_raw_L2.shape[0]>1:
				Repeat_raw_L2_sort = remove_redundancy.clean_df(Repeat_raw_L2)
				Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2_sort,ignore_index = True)
			elif Repeat_raw_L2.shape[0]==1:
				Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2,ignore_index = True)

	if Repeat_df_L2.shape[0]>1:
		Repeat_df_L2=remove_redundancy.clean_df(Repeat_df_L2)
	return Repeat_df_L2

def sliding_window_micro(window_size,V_start,V_length,b,based3):
	Repeat_df_L3 = pd.DataFrame(columns=["Motifs", "Period","RepeatS","Size","Start","End","Position_O","AlignmentS","GapS"])
	for j in range(V_start-window_size+1,V_start+V_length): 
		V_pattern = b[j:j+window_size]
		Repeat_raw_L3 = repeatAlign.count_microsatellite(j,V_pattern,window_size,b,V_start,based3)
		if Repeat_raw_L3 != None:
			Repeat_df_L3.loc[j]=Repeat_raw_L3
	return Repeat_df_L3

def sliding_window_mini(window_size,V_start,V_length,b,alignAS3,gapAS3,based3):
	Repeat_df_L3 = pd.DataFrame(columns=["Motifs", "Period","RepeatS","Size","Start","End","Position_O","AlignmentS","GapS"])
	for j in range(V_start-window_size+1,V_start+V_length):
		V_pattern = b[j:j+window_size]
		Repeat_raw_L3 = repeatAlign.count_minisatellite(j,V_pattern,window_size,b,V_start,alignAS3,gapAS3,based3)
		if Repeat_raw_L3 != None:
			Repeat_df_L3.loc[j]=Repeat_raw_L3
	return Repeat_df_L3

def search_repeats(df,ref_file,min_size0,max_size0,exact_size0,alignAS0,gapAS0,based0):
	Repeat_df_L0S = pd.DataFrame(columns=["Chromosome","Position","REF","ALT","Genotype","Motifs","Period","Size","Start","End","Repeat_Score","Alignment_Score","Gap_Score","Repeat_GC%"])
	
	chr_name=""
	for i in range(0,df.shape[0]):
		chromosome =  df.loc[i,"Chromosome"]
		position = df.loc[i,"Position"]
		REF = df.loc[i,"REF"]
		ALT = df.loc[i,"ALT"]
		Genotype = df.loc[i,"Genotype"]

		if chr_name!=chromosome:
			chr_name=chromosome
			Ref=read_reference.read_reference_sequence(ref_file,chr_name)

		Repeat_raw_L0 = catching_variant(position,REF,ALT,Ref,min_size0,max_size0,exact_size0,alignAS0,gapAS0,based0) 

		Repeat_df_L0 = pd.DataFrame(columns=["Chromosome","Position","REF","ALT","Genotype","Motifs","Period","Size","Start","End","Repeat_Score","Alignment_Score","Gap_Score","Repeat_GC%"])
		if Repeat_raw_L0.empty == False:
			for j in range(0,Repeat_raw_L0.shape[0]):
				GC_Percent = Local_GC.GC_content(Ref[Repeat_raw_L0.iloc[j]["Start"]-1:Repeat_raw_L0.iloc[j]["End"]])
				item=[chromosome,position,REF,ALT,Genotype,str(Repeat_raw_L0.iloc[j]["Motifs"]),Repeat_raw_L0.iloc[j]["Period"],Repeat_raw_L0.iloc[j]["Size"],Repeat_raw_L0.iloc[j]["Start"],Repeat_raw_L0.iloc[j]["End"],round(Repeat_raw_L0.iloc[j]["RepeatS"],2),round(Repeat_raw_L0.iloc[j]["AlignmentS"],2),round(Repeat_raw_L0.iloc[j]["GapS"],2),GC_Percent]
				Repeat_df_L0.loc[j]=item

			Repeat_df_L0S = Repeat_df_L0S.append(Repeat_df_L0,ignore_index = True)
	
		else:
			continue
		
	return Repeat_df_L0S



































