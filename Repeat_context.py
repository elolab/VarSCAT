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

def catching_variant(position_l1, REF_l1, ALT_l1, Ref1,min_size1,max_size1,min_time1,min_per1,alignAS1,based1,Mscore1,MIscore1,OGscore1,Rcheck1,S_except1):
	if REF_l1[0]==ALT_l1[0]:
		if len(REF_l1) > len(ALT_l1):
			variant_start = int(position_l1)-based1
			variant_length = len(REF_l1)+1
		elif len(REF_l1) < len(ALT_l1):
			variant_start = int(position_l1)-based1
			variant_length = 2 
	elif REF_l1[0]!=ALT_l1[0]:
		if len(REF_l1) > len(ALT_l1):
			variant_start = int(position_l1)-based1-1
			variant_length = len(REF_l1)+2
		elif len(REF_l1) < len(ALT_l1):
			variant_start = int(position_l1)-based1-1
			variant_length = len(REF_l1)+2 
		elif len(REF_l1) == len(ALT_l1):
			variant_start = int(position_l1)-based1-1
			variant_length = len(REF_l1)+2 

	Repeat_raw_L1 = sizing_window(variant_start,variant_length,Ref1,min_size1,max_size1,min_time1,min_per1,alignAS1,based1,Mscore1,MIscore1,OGscore1,Rcheck1,S_except1)

	return Repeat_raw_L1

def sizing_window(start,length,Ref2,min_size2,max_size2,min_time2,min_per2,alignAS2,based2,Mscore2,MIscore2,OGscore2,Rcheck2,S_except2):
	Repeat_df_L2 = pd.DataFrame(columns=["Motifs", "Period","Size","Start","End","RepeatS","AlignmentS","Match","Mismatch","Gap"])
	max_RS2=0
	for i in range(min_size2,max_size2+1): 
		mini_wiondow = sliding_window_mini(i,start,length,Ref2,alignAS2,based2,min_time2,min_per2,Mscore2,MIscore2,OGscore2,Rcheck2,S_except2)
		Repeat_raw_L2 = mini_wiondow[0]
		if Repeat_raw_L2.shape[0]>1:
			Repeat_raw_L2_sort = remove_redundancy.clean_df(Repeat_raw_L2,Ref2,based2)
			Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2_sort,ignore_index = True)
		elif Repeat_raw_L2.shape[0]==1:
			Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2,ignore_index = True)
		elif Repeat_raw_L2.empty == True:
			if mini_wiondow[1]>max_RS2:
				max_RS2=mini_wiondow[1]
				trim_list2=mini_wiondow[2]

	if Repeat_df_L2.shape[0]>1:
		Repeat_df_L2=remove_redundancy.clean_df(Repeat_df_L2,Ref2,based2)
	elif Repeat_df_L2.empty == True and max_RS2!=0:
		repeat_trimmed=repeatAlign.trim_repeat(trim_list2,min_time2,min_per2,Mscore2,MIscore2,OGscore2,alignAS2,based2)
		if repeat_trimmed.shape[0]==1:
			Repeat_df_L2=Repeat_df_L2.append(repeat_trimmed,ignore_index = True)
		
	return Repeat_df_L2

def sliding_window_mini(window_size,V_start,V_length,b,alignAS3,based3,min_time3,min_per3,Mscore3,MIscore3,OGscore3,Rcheck3,S_except3):
	Repeat_df_L3 = pd.DataFrame(columns=["Motifs", "Period","Size","Start","End","RepeatS","AlignmentS","Match","Mismatch","Gap"])
	max_RS3=0
	trim_list3=""
	for j in range(V_start-window_size+1,V_start+V_length):
		V_pattern = b[j:j+window_size]
		if "N" not in V_pattern:
			Repeat_raw_L3 = repeatAlign.count_minisatellite(j,V_pattern,window_size,b,based3,min_time3,Mscore3,MIscore3,OGscore3,Rcheck3,S_except3)
		else:
			continue
		if Repeat_raw_L3 != None:
			if Repeat_raw_L3[6] >= alignAS3 and Repeat_raw_L3[7] >= min_per3:
				Repeat_df_L3.loc[j]=Repeat_raw_L3[0:10]
			else:
				if Repeat_raw_L3[7]>max_RS3:
					max_RS3=Repeat_raw_L3[6]
					trim_list3=Repeat_raw_L3

	return Repeat_df_L3,max_RS3,trim_list3

def search_repeats(df,ref_file,min_size0,max_size0,min_time0,min_per0,alignAS0,based0,Mscore0,MIscore0,OGscore0,Rcheck0,S_except0):
	chr_name=""
	Motifs_list=[]
	Period_list=[]
	Size_list=[]
	Start_list=[]
	End_list=[]
	Repeat_Score_list=[]
	Alignment_Score_list=[]
	Match_list=[]
	Mismatch_list=[]
	Gap_list=[]
	Repeat_GC_list=[]
	for i in range(0,df.shape[0]):
		chromosome =  df.loc[i,"Chromosome"]
		position = df.loc[i,"Position"]
		REF = df.loc[i,"REF"]
		ALT = df.loc[i,"ALT"]

		Motifs_string=""
		Period_string=""
		Size_string=""
		Start_string=""
		End_string=""
		Repeat_Score_string=""
		Alignment_Score_string=""
		Match_string=""
		Mismatch_string=""
		Gap_string=""
		Repeat_GC_string=""

		if ("N" not in REF) and ("N" not in ALT): 

			if chr_name!=chromosome:
				chr_name=chromosome
				Ref=read_reference.read_reference_sequence(ref_file,chr_name)

			Repeat_raw_L0 = catching_variant(position,REF,ALT,Ref,min_size0,max_size0,min_time0,min_per0,alignAS0,based0,Mscore0,MIscore0,OGscore0,Rcheck0,S_except0) 

			if Repeat_raw_L0.empty == False:
				for j in range(0,Repeat_raw_L0.shape[0]):
					if Repeat_raw_L0.iloc[j]["Motifs"].count("A")==Repeat_raw_L0.iloc[j]["Size"] or Repeat_raw_L0.iloc[j]["Motifs"].count("T")==Repeat_raw_L0.iloc[j]["Size"] or Repeat_raw_L0.iloc[j]["Motifs"].count("G")==Repeat_raw_L0.iloc[j]["Size"] or Repeat_raw_L0.iloc[j]["Motifs"].count("C")==Repeat_raw_L0.iloc[j]["Size"] and Repeat_raw_L0.iloc[j]["Size"]!=1:
						Motifs_string=Motifs_string+str(Repeat_raw_L0.iloc[j]["Motifs"][0])
						Period_string=Period_string+str(int(Repeat_raw_L0.iloc[j]["Period"]*Repeat_raw_L0.iloc[j]["Size"]))
						Size_string=Size_string+str(1)
						Start_string=Start_string+str(int(Repeat_raw_L0.iloc[j]["Start"]))
						End_string=End_string+str(int(Repeat_raw_L0.iloc[j]["End"]))
						Repeat_Score_string=Repeat_Score_string+str(round((Repeat_raw_L0.iloc[j]["AlignmentS"]/(Repeat_raw_L0.iloc[j]["End"]-Repeat_raw_L0.iloc[j]["Start"]+1))*(Repeat_raw_L0.iloc[j]["Period"]*Repeat_raw_L0.iloc[j]["Size"]),2))
						Alignment_Score_string=Alignment_Score_string+str(Repeat_raw_L0.iloc[j]["AlignmentS"])
						Match_string=Match_string+str(round(Repeat_raw_L0.iloc[j]["Match"]*100,2))
						Mismatch_string=Mismatch_string+str(round(Repeat_raw_L0.iloc[j]["Mismatch"]*100,2))
						Gap_string=Gap_string+str(round(Repeat_raw_L0.iloc[j]["Gap"]*100,2))
						GC_Percent = Local_GC.GC_content(Ref[int(Repeat_raw_L0.iloc[j]["Start"]-based0):int(Repeat_raw_L0.iloc[j]["End"]-based0+1)].upper())
						Repeat_GC_string=Repeat_GC_string+str(GC_Percent)
					else:
						Motifs_string=Motifs_string+str(Repeat_raw_L0.iloc[j]["Motifs"])
						Period_string=Period_string+str(int(Repeat_raw_L0.iloc[j]["Period"]))
						Size_string=Size_string+str(int(Repeat_raw_L0.iloc[j]["Size"]))
						Start_string=Start_string+str(int(Repeat_raw_L0.iloc[j]["Start"]))
						End_string=End_string+str(int(Repeat_raw_L0.iloc[j]["End"]))
						Repeat_Score_string=Repeat_Score_string+str(round(Repeat_raw_L0.iloc[j]["RepeatS"],2))
						Alignment_Score_string=Alignment_Score_string+str(Repeat_raw_L0.iloc[j]["AlignmentS"])
						Match_string=Match_string+str(round(Repeat_raw_L0.iloc[j]["Match"]*100,2))
						Mismatch_string=Mismatch_string+str(round(Repeat_raw_L0.iloc[j]["Mismatch"]*100,2))
						Gap_string=Gap_string+str(round(Repeat_raw_L0.iloc[j]["Gap"]*100,2))
						GC_Percent = Local_GC.GC_content(Ref[int(Repeat_raw_L0.iloc[j]["Start"]-based0):int(Repeat_raw_L0.iloc[j]["End"]-based0+1)].upper())
						Repeat_GC_string=Repeat_GC_string+str(GC_Percent)

					if j<(Repeat_raw_L0.shape[0]-1):
						Motifs_string=Motifs_string+","
						Period_string=Period_string+","
						Size_string=Size_string+","
						Start_string=Start_string+","
						End_string=End_string+","
						Repeat_Score_string=Repeat_Score_string+","
						Alignment_Score_string=Alignment_Score_string+","
						Match_string=Match_string+","
						Mismatch_string=Mismatch_string+","
						Gap_string=Gap_string+","
						GC_Percent = Local_GC.GC_content(Ref[int(Repeat_raw_L0.iloc[j]["Start"]-based0):int(Repeat_raw_L0.iloc[j]["End"]-based0+1)].upper())
						Repeat_GC_string=Repeat_GC_string+","
			
				Motifs_list.append(Motifs_string)
				Period_list.append(Period_string)
				Size_list.append(Size_string)
				Start_list.append(Start_string)
				End_list.append(End_string)
				Repeat_Score_list.append(Repeat_Score_string)
				Alignment_Score_list.append(Alignment_Score_string)
				Match_list.append(Match_string)
				Mismatch_list.append(Mismatch_string)
				Gap_list.append(Gap_string)
				Repeat_GC_list.append(Repeat_GC_string)
	
			else:
				Motifs_list.append(Motifs_string)
				Period_list.append(Period_string)
				Size_list.append(Size_string)
				Start_list.append(Start_string)
				End_list.append(End_string)
				Repeat_Score_list.append(Repeat_Score_string)
				Alignment_Score_list.append(Alignment_Score_string)
				Match_list.append(Match_string)
				Mismatch_list.append(Mismatch_string)
				Gap_list.append(Gap_string)
				Repeat_GC_list.append(Repeat_GC_string)

		else:
			Motifs_list.append(Motifs_string)
			Period_list.append(Period_string)
			Size_list.append(Size_string)
			Start_list.append(Start_string)
			End_list.append(End_string)
			Repeat_Score_list.append(Repeat_Score_string)
			Alignment_Score_list.append(Alignment_Score_string)
			Match_list.append(Match_string)
			Mismatch_list.append(Mismatch_string)
			Gap_list.append(Gap_string)
			Repeat_GC_list.append(Repeat_GC_string)
	
	return Motifs_list,Period_list,Size_list,Start_list,End_list,Repeat_Score_list,Alignment_Score_list,Match_list,Mismatch_list,Gap_list,Repeat_GC_list

