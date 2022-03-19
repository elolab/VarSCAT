from Bio.Seq import Seq
from Bio import Align
from Bio import motifs
import pandas as pd
 
def alignment_pattern(PS,SZ,RF,PR,DIR,SZ_E): 
	aligner = Align.PairwiseAligner()
	aligner.mode = "global"
	aligner.match_score = 5
	aligner.mismatch_score = -5
	aligner.open_gap_score = -1000.0 

	if type(SZ_E)!=int:
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
		score_R = aligner.score(PR, potential_R_N)
		
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
	num_mismatch_R = [0]
	num_gap_L = []
	num_mismatch_L = []
	
	E_L=[repeat_end]
	S_L=[repeat_start]
	
	repeat_pass =5*repeat_check-5*(1-repeat_check)
	Align_inner = alignment_pattern(pattern_start,search_size,ref,potential_R,1,size_except)
	while float(Align_inner[0])/len(potential_R)>=repeat_pass:  
		if "N" not in Align_inner[2]:
			memory_pattern_list.append(Align_inner[2])
		
			num_gap_R.append(Align_inner[1])
			mismatch_num= ((5*search_size)-Align_inner[0])/10
			num_mismatch_R.append(mismatch_num)
			
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
	while float(Align_inner[0])/len(potential_R) >=repeat_pass:
		if "N" not in Align_inner[2]:
			memory_pattern_list.append(Align_inner[2])

			num_gap_L.append(Align_inner[1])
			mismatch_num= ((5*search_size)-Align_inner[0])/10
			num_mismatch_L.append(mismatch_num)			
			
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
		num_mismatch=num_mismatch_R+num_mismatch_L

		aligner = Align.PairwiseAligner()
		aligner.mode = "global"
		aligner.match_score = 5
		aligner.mismatch_score = -5
		aligner.open_gap_score = -1000.0
		score_R = aligner.score(potential_R0, potential_R)
		mismatch_num= ((5*search_size)-score_R)/10
		num_mismatch.insert(0,mismatch_num)
		for i in range(1,4):
			score_R = aligner.score(memory_pattern_list[i], potential_R)
			mismatch_num= ((5*search_size)-score_R)/10
			num_mismatch.insert(1,mismatch_num)
		del num_mismatch[4:8]
		num_mismatch=num_mismatch[len(num_mismatch_R):][::-1]+num_mismatch[0:len(num_mismatch_R)]

		sum_mismatch=sum(num_mismatch)
		sum_gap=sum(num_gap)
		T_MS=float(sum_mismatch)/repeat_length
		T_G=float(sum_gap)/repeat_length
		T_M=float(repeat_length-sum_mismatch-sum_gap)/repeat_length

		AS=(repeat_length-sum_mismatch-sum_gap)*match_score+sum_mismatch*mismatch_score+sum_gap*gap_score
		repeat_score = (float(AS)/(repeat_length))*repeat_time		

		return potential_R.upper(),repeat_time,len(potential_R.upper()),repeat_start+based,repeat_end+based,repeat_score,AS,T_M,T_MS,T_G,S_L,E_L,memory_pattern_list,num_mismatch,num_gap

def trim_repeat(RL,MT,MP,M_score,MIS_score,G_score,ASC,BS):
	R_DF = pd.DataFrame(columns=["Motifs", "Period","Size","Start","End","RepeatS","AlignmentS","Match","Mismatch","Gap"])
	S_L=RL[10]
	E_L=RL[11]
	S_L1=S_L[0] 
	E_L1=E_L[0]
	MM_R=RL[13][len(S_L)-1:] 
	MM_L=RL[13][0:len(S_L)-1]
	MM_L.append(MM_R[0])
	GP_R=RL[14][len(S_L)-1:] 
	GP_L=RL[14][0:len(S_L)-1]
	RT=RL[1]
	P_R=RL[12][0:len(E_L)]
	P_L=RL[12][len(E_L):]
	while len(GP_R)>0 and RT>=MT:
		R_len=E_L[-1]-S_L1+1
		MA=float(R_len-sum(MM_R)-sum(GP_R))/R_len	
		if MA < MP: 
			del MM_R[-1]		
			RT=RT-1	
			del GP_R[-1]		
			P_R=P_R[:-1]
			del E_L[-1]	
		else:
			break
	
	while len(GP_L)>0 and RT>=MT:
		R_len=E_L1-S_L[-1]+1
		MA=float(R_len-sum(MM_L)-sum(GP_L))/R_len	
		if MA < MP:
			del MM_L[0]
			RT=RT-1
			del GP_L[0]
			P_L=P_L[:-1]
			del S_L[-1]
		else:
			break

	if len(MM_L)>0:
		del MM_L[-1]
	
	R_len=E_L[-1]-S_L[-1]+1
	
	if RT >= MT:
		sum_mismatch=sum(MM_R+MM_L)
		sum_gap=sum(GP_R+GP_L)
		MA=float(R_len-sum_mismatch-sum_gap)/R_len	
		AS=(R_len-sum_mismatch-sum_gap)*M_score+sum_mismatch*MIS_score+sum_gap*G_score
		if AS >= ASC and MA >= MP:	 
			MM=float(sum_mismatch)/R_len
			GP=float(sum_gap)/R_len	
			RS = (float(AS)/R_len)*RT
			m = motifs.create(P_R+P_L)
			m1 = m.consensus
			R_DF.loc[0]=(m1.upper(),RT,RL[2],S_L[-1]+BS,E_L[-1]+BS,RS,AS,MA,MM,GP)
	
	return R_DF

