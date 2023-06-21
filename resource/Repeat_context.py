from .repeatAlign import count_minisatellite, trim_repeat
from .remove_redundancy import clean_df
from .Local_GC import GC_content
from .read_reference import read_reference_sequence
from .copy_change import check_copy

def sizing_window(start,end,R_string2,A_string2,Ref2,min_size2,max_size2,min_time2,min_per2,alignAS2,based2,Mscore2,MIscore2,OGscore2,Rcheck2,S_except2):
	# variant site should be equivalent regions
	if R_string2[0]==A_string2[0] and len(R_string2) < len(A_string2):
		start = int(start)-based2
	else:
		start = int(start)-based2-1
	end = int(end)-based2+1
	
	Repeat_df_L2 = []
	trim_list2= []
	for i in range(min_size2,max_size2+1): 
		mini_wiondow = sliding_window_mini(i,start,end,Ref2,alignAS2,based2,min_time2,min_per2,Mscore2,MIscore2,OGscore2,Rcheck2,S_except2)
		Repeat_raw_L2 = mini_wiondow[0]
		if len(Repeat_raw_L2)>1:
			Repeat_raw_L2_sort = clean_df(Repeat_raw_L2,Ref2,based2)
			Repeat_df_L2=Repeat_df_L2+Repeat_raw_L2_sort
		elif len(Repeat_raw_L2)==1:
			Repeat_df_L2=Repeat_df_L2+Repeat_raw_L2
		elif len(Repeat_raw_L2) == 0 and mini_wiondow[1]!=0:
			trim_list2.append(mini_wiondow[2])

	if len(Repeat_df_L2)>1:
		Repeat_df_L2=clean_df(Repeat_df_L2,Ref2,based2)
	elif len(Repeat_df_L2) == 0 and len(trim_list2)!=0:
		for candidate in trim_list2:
			repeat_trimmed=trim_repeat(candidate,min_time2,min_per2,Mscore2,MIscore2,OGscore2,alignAS2,based2)
			if len(repeat_trimmed)!=0:
				Repeat_df_L2.append(repeat_trimmed)
		if len(Repeat_df_L2)>1:
			Repeat_df_L2=clean_df(Repeat_df_L2,Ref2,based2)
		
	return Repeat_df_L2

def sliding_window_mini(window_size,V_start,V_end,b,alignAS3,based3,min_time3,min_per3,Mscore3,MIscore3,OGscore3,Rcheck3,S_except3):
	Repeat_df_L3 = []
	max_RS3=0
	trim_list3=""
	V_memo=[]
	for j in range(V_start-window_size+1,V_end+1):
		V_pattern = b[j:j+window_size]
		
		# skip if the same pattern have been searched
		if j<=V_start:# This means pattern is not in permutated region anymore
			V_memo.append(V_pattern)
		elif j>V_start and j<=V_end-window_size: # # This means pattern is in permutated region anymore, j>V_end-window_size means not in premutated region
			# check if ther premutated pattern have been searched	
			if V_pattern == V_memo[len(V_memo)-window_size]:
				continue
			else:
				V_memo.append(V_pattern)
			
			
		if "N" not in V_pattern:
			Repeat_raw_L3 = count_minisatellite(j,V_pattern,window_size,b,based3,min_time3,Mscore3,MIscore3,OGscore3,Rcheck3,S_except3)
		else:
			continue
		if Repeat_raw_L3 != None:
			if Repeat_raw_L3[6] >= alignAS3 and Repeat_raw_L3[7] >= min_per3:
				Repeat_df_L3.append(Repeat_raw_L3[0:10])
			else:
				if Repeat_raw_L3[6]>max_RS3:
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
	copynumber_change=[]
	for i in range(0,len(df)):
		chromosome =  df[i][0]
		REF = df[i][2]
		ALT = df[i][3]
		left_align = df[i][4]
		right_edge = df[i][6]

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
		copynumber_change_string=""

		if ("N" not in REF) and ("N" not in ALT): 

			if chr_name!=chromosome:
				chr_name=chromosome
				Ref=read_reference_sequence(ref_file,chr_name)

			Repeat_raw_L0 = sizing_window(left_align,right_edge,REF,ALT,Ref,min_size0,max_size0,min_time0,min_per0,alignAS0,based0,Mscore0,MIscore0,OGscore0,Rcheck0,S_except0) 

			if len(Repeat_raw_L0) != 0:
				for j in range(0,len(Repeat_raw_L0)):
					if Repeat_raw_L0[j][0].count("A")==Repeat_raw_L0[j][2] or Repeat_raw_L0[j][0].count("T")==Repeat_raw_L0[j][2] or Repeat_raw_L0[j][0].count("G")==Repeat_raw_L0[j][2] or Repeat_raw_L0[j][0].count("C")==Repeat_raw_L0[j][2] and Repeat_raw_L0[j][2]!=1:
						Motifs_string=Motifs_string+str(Repeat_raw_L0[j][0][0])
						Period_string=Period_string+str(int(Repeat_raw_L0[j][1]*Repeat_raw_L0[j][2]))
						Size_string=Size_string+str(1)
						Start_string=Start_string+str(int(Repeat_raw_L0[j][3]))
						End_string=End_string+str(int(Repeat_raw_L0[j][4]))
						Repeat_Score_string=Repeat_Score_string+str(round((Repeat_raw_L0[j][6]/(Repeat_raw_L0[j][4]-Repeat_raw_L0[j][3]+1))*(Repeat_raw_L0[j][1]*Repeat_raw_L0[j][2]),2))
						Alignment_Score_string=Alignment_Score_string+str(Repeat_raw_L0[j][6])
						Match_string=Match_string+str(round(Repeat_raw_L0[j][7]*100,2))
						Mismatch_string=Mismatch_string+str(round(Repeat_raw_L0[j][8]*100,2))
						Gap_string=Gap_string+str(round(Repeat_raw_L0[j][9]*100,2))
						GC_Percent = GC_content(Ref[int(Repeat_raw_L0[j][3]-based0):int(Repeat_raw_L0[j][4]-based0+1)].upper())
						Repeat_GC_string=Repeat_GC_string+str(GC_Percent)
						cp_change=check_copy(REF,ALT,int(left_align),int(right_edge),int(Repeat_raw_L0[j][3]),int(Repeat_raw_L0[j][4]),str(Repeat_raw_L0[j][0][0]),int(Repeat_raw_L0[j][1]*Repeat_raw_L0[j][2]),Repeat_raw_L0[j][7],Rcheck0,S_except0)
						copynumber_change_string=copynumber_change_string+str(cp_change)
					else:
						Motifs_string=Motifs_string+str(Repeat_raw_L0[j][0])
						Period_string=Period_string+str(int(Repeat_raw_L0[j][1]))
						Size_string=Size_string+str(int(Repeat_raw_L0[j][2]))
						Start_string=Start_string+str(int(Repeat_raw_L0[j][3]))
						End_string=End_string+str(int(Repeat_raw_L0[j][4]))
						Repeat_Score_string=Repeat_Score_string+str(round(Repeat_raw_L0[j][5],2))
						Alignment_Score_string=Alignment_Score_string+str(Repeat_raw_L0[j][6])
						Match_string=Match_string+str(round(Repeat_raw_L0[j][7]*100,2))
						Mismatch_string=Mismatch_string+str(round(Repeat_raw_L0[j][8]*100,2))
						Gap_string=Gap_string+str(round(Repeat_raw_L0[j][9]*100,2))
						GC_Percent = GC_content(Ref[int(Repeat_raw_L0[j][3]-based0):int(Repeat_raw_L0[j][4]-based0+1)].upper())
						Repeat_GC_string=Repeat_GC_string+str(GC_Percent)
						cp_change=check_copy(REF,ALT,int(left_align),int(right_edge),int(Repeat_raw_L0[j][3]),int(Repeat_raw_L0[j][4]),str(Repeat_raw_L0[j][0]),int(Repeat_raw_L0[j][1]),Repeat_raw_L0[j][7],Rcheck0,S_except0)
						copynumber_change_string=copynumber_change_string+str(cp_change)

					if j<(len(Repeat_raw_L0)-1):
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
						#GC_Percent = GC_content(Ref[int(Repeat_raw_L0[j][3]-based0):int(Repeat_raw_L0[j][4]-based0+1)].upper())
						Repeat_GC_string=Repeat_GC_string+","
						copynumber_change_string=copynumber_change_string+","
			
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
				copynumber_change.append(copynumber_change_string)
	
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
				copynumber_change.append(copynumber_change_string)

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
			copynumber_change.append(copynumber_change_string)
	
	return Motifs_list,Period_list,Size_list,Start_list,End_list,Repeat_Score_list,Alignment_Score_list,Match_list,Mismatch_list,Gap_list,Repeat_GC_list,copynumber_change

