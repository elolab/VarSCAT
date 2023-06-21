from Bio.Seq import Seq
from Bio import motifs

def diff_letters(a,b):
	common = sum(a[i]==b[i] for i in range(len(a)))
	return common

def alignment_pattern(PS,SZ,RF,PR,DIR,SZ_E): 
	if SZ_E == -1:
		SZ_assign=SZ+1
	else:
		SZ_assign=SZ_E+1

	list_loc = []
	list_score = []
	list_seq = []
	for i in range(0,SZ_assign):
		if DIR == 1:
			potential_R_N = Seq(RF[PS+SZ+i:PS+SZ+SZ+i].upper())
		elif DIR == -1:
			potential_R_N = Seq(RF[PS-SZ-i:PS-i].upper())
		score_R = diff_letters(PR, potential_R_N)
		
		list_loc.append(i)
		list_score.append(score_R)
		list_seq.append(potential_R_N)

	MS_score = max(list_score)
	MS_index = list_score.index(MS_score)
	
	MS_distance = list_loc[MS_index]
	MS_seq = list_seq[MS_index]

	return MS_score, MS_distance, MS_seq

def count_minisatellite(pattern_start,potential_R,search_size,ref,based,mperiod,match_score,mismatch_score,gap_score,repeat_check,size_except):
	repeat_time = 1

	repeat_start = pattern_start
	repeat_end = pattern_start + search_size -1

	repeat_start1 = repeat_start
	repeat_end1 = repeat_end

	pattern_start1 = pattern_start

	potential_R = Seq(potential_R.upper())
	potential_R0 = potential_R 
	memory_pattern_list = [potential_R] 

	num_gap_R = []
	num_gap_L = []
	
	E_L=[repeat_end]
	S_L=[repeat_start]
	
	Align_inner = alignment_pattern(pattern_start,search_size,ref,potential_R,1,size_except)
	while float(Align_inner[0])/len(potential_R)>=repeat_check:  
		if "N" not in Align_inner[2]:
			memory_pattern_list.append(Align_inner[2])
		
			num_gap_R.append(Align_inner[1])
			
			if len(memory_pattern_list) > 4: 
				m = motifs.create(memory_pattern_list)
				potential_R = m.consensus
			else:
				potential_R = potential_R
			pattern_start = pattern_start+search_size + Align_inner[1] 
			repeat_time = repeat_time +1
			repeat_end = repeat_end + search_size + Align_inner[1]
			E_L.append(repeat_end)
		
			Align_inner = alignment_pattern(pattern_start,search_size,ref,potential_R,1,size_except)				
		else:
			break
 
	Align_inner = alignment_pattern(pattern_start1,search_size,ref,potential_R,-1,size_except)  
	while float(Align_inner[0])/len(potential_R) >=repeat_check:
		if "N" not in Align_inner[2]:
			memory_pattern_list.append(Align_inner[2])

			num_gap_L.append(Align_inner[1])		
			
			if len(memory_pattern_list) > 4:
				m = motifs.create(memory_pattern_list)
				potential_R = m.consensus
			else:
				potential_R = potential_R
			pattern_start1 = pattern_start1-search_size-Align_inner[1]
			repeat_time = repeat_time +1
			repeat_start = repeat_start-search_size-Align_inner[1]
			S_L.append(repeat_start)

			Align_inner = alignment_pattern(pattern_start1,search_size,ref,potential_R,-1,size_except)
		else:
			break

	repeat_length = repeat_end-repeat_start+1

	if (repeat_start != repeat_start1 or repeat_end != repeat_end1) and repeat_time >=mperiod:
		num_gap=num_gap_L[::-1]+num_gap_R
		memory_pattern_list=memory_pattern_list[len(num_gap_R)+1:][::-1]+memory_pattern_list[0:len(num_gap_R)+1]
		num_mismatch=[]
		for i in range(0,len(memory_pattern_list)):
			score_R = diff_letters(memory_pattern_list[i], potential_R)
			mismatch_num= search_size-score_R
			num_mismatch.append(mismatch_num)
		
		sum_mismatch=sum(num_mismatch)
		sum_gap=sum(num_gap)
		T_MS=float(sum_mismatch)/repeat_length
		T_G=float(sum_gap)/repeat_length
		T_M=float(repeat_length-sum_mismatch-sum_gap)/repeat_length

		AS=(repeat_length-sum_mismatch-sum_gap)*match_score+sum_mismatch*mismatch_score+sum_gap*gap_score
		repeat_score = (float(AS)/(repeat_length))*repeat_time		

		return potential_R.upper(),repeat_time,len(potential_R.upper()),repeat_start+based,repeat_end+based,repeat_score,AS,T_M,T_MS,T_G,S_L,E_L,memory_pattern_list,num_mismatch,num_gap

def trim_repeat(RL,MT,MP,M_score,MIS_score,G_score,ASC,BS):
	R_DF = []
	S_L=RL[10][::-1]
	E_L=RL[11]
	L_marker=["S"]*len(S_L)+["E"]*len(E_L)
	L_L=S_L+E_L
	MM_L=RL[13]
	GP_L=RL[14]
	P_L=RL[12]
	RT=RL[1]
	
	while len(GP_L)>0 and RT>=MT:
		if L_marker[-1]=="E" and L_marker[0]=="S":
			R_len=L_L[-1]-L_L[0]+1
		elif L_marker[-1]=="E" and L_marker[0]=="E":
			R_len=L_L[-1]-L_L[0]
		elif L_marker[-1]=="S" and L_marker[0]=="S":
			R_len=L_L[-1]-L_L[0]
		
		MA=float(R_len-sum(MM_L)-sum(GP_L))/R_len
		if MA >= MP:
			sum_mismatch=sum(MM_L)
			sum_gap=sum(GP_L)
			AS=(R_len-sum_mismatch-sum_gap)*M_score+sum_mismatch*MIS_score+sum_gap*G_score
			if AS >=ASC:
				MM=float(sum_mismatch)/R_len
				GP=float(sum_gap)/R_len	
				RS = (float(AS)/R_len)*RT
				m = motifs.create(P_L)
				m1 = m.consensus
				if L_marker[0] == "E":
					 L_L[0]=L_L[0]+1
				if L_marker[-1] == "S":
					 L_L[-1]=L_L[-1]-1
				R_DF=[m1.upper(),RT,RL[2],L_L[0]+BS,L_L[-1]+BS,RS,AS,MA,MM,GP]
				
				break
		
		next=0
		while next<=(len(GP_L)/2):
			penalty_L=MM_L[next]*MIS_score+GP_L[next]*G_score
			penalty_R=MM_L[-next-1]*MIS_score+GP_L[-next-1]*G_score
			if penalty_L > penalty_R:
				del MM_L[-1]		
				RT=RT-1	
				del GP_L[-1]		
				P_L=P_L[:-1]
				del L_L[-1]
				del L_marker[-1]
				break
			elif penalty_L < penalty_R:
				del MM_L[0]
				RT=RT-1
				del GP_L[0]
				P_L=P_L[1:]
				del L_L[0]
				del L_marker[0]
				break
			next=next+1
		if next>(len(GP_L)/2):
			del MM_L[-1]		
			RT=RT-1	
			del GP_L[-1]		
			P_L=P_L[:-1]
			del L_L[-1]
			del L_marker[-1]
	
	return R_DF
	
