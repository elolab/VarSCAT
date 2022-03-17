from Bio.Seq import Seq
from Bio import Align
from Bio import motifs


def alignment_pattern(PS,SZ,RF,PR,DIR): 
	aligner = Align.PairwiseAligner()
	aligner.mode = "global"
	aligner.match_score = 5.0
	aligner.mismatch_score = -4.0
	aligner.open_gap_score = -50.0 

	list_loc = []
	list_score = []
	list_seq = []
	for i in range(0,SZ):
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
	if MS_distance==0:
		MS_distance_score = float(MS_distance)
	else:
		MS_distance_score = float(10+(MS_distance-1)*1)
	MS_seq = list_seq[MS_index]

	return MS_score, MS_distance, MS_seq, MS_distance_score



def count_minisatellite(pattern_start, potential_R, search_size, ref, variant_start,alignAS,gapAS,based):
	
	repeat_time = 1
	repeat_start = pattern_start
	repeat_end = pattern_start + search_size -1

	repeat_start1 = repeat_start
	repeat_end1 = repeat_end

	pattern_start1 = pattern_start

	potential_R = Seq(potential_R.upper())
	potential_R0 = potential_R 
	memory_pattern_list_R = [] 
	memory_pattern_list_R.append(potential_R)

	memory_distance_list = []
	memory_score_list = []

	Align_inner = alignment_pattern(pattern_start,search_size,ref,potential_R,1)  
	while float(Align_inner[0])/len(potential_R) >=2.75:  
		if "N" not in Align_inner[2]:
			memory_pattern_list_R.append(Align_inner[2])
			memory_distance_list.append(Align_inner[3]) 
			memory_score_list.append(Align_inner[0])
			if len(memory_pattern_list_R) > 4: 
				m = motifs.create(memory_pattern_list_R)
				potential_R = m.consensus
			else:
				potential_R = potential_R
			pattern_start = pattern_start+search_size + Align_inner[1] 
			repeat_time = repeat_time +1
			repeat_end = repeat_end + search_size + Align_inner[1]
			Align_inner = alignment_pattern(pattern_start,search_size,ref,potential_R,1)
		else:
			break

	Align_inner = alignment_pattern(pattern_start1,search_size,ref,potential_R,-1)  
	while float(Align_inner[0])/len(potential_R) >=2.75:
		if "N" not in Align_inner[2]:
			memory_pattern_list_R.append(Align_inner[2])
			memory_distance_list.append(Align_inner[3])
			memory_score_list.append(Align_inner[0])
			if len(memory_pattern_list_R) > 4:
				m = motifs.create(memory_pattern_list_R)
				potential_R = m.consensus
			else:
				potential_R = potential_R
			pattern_start1 = pattern_start1 - search_size -  Align_inner[1]
			repeat_time = repeat_time +1
			repeat_start = repeat_start -search_size -  Align_inner[1]
			Align_inner = alignment_pattern(pattern_start1,search_size,ref,potential_R,-1)
		else:
			break

	if repeat_start != repeat_start1 or repeat_end != repeat_end1:
		if repeat_time >=5:
			aligner = Align.PairwiseAligner()
			aligner.mode = "global"
			aligner.match_score = 5.0
			aligner.mismatch_score = -4.0
			aligner.open_gap_score = -50
			score_R = aligner.score(potential_R0, potential_R)
			memory_score_list.append(score_R)
			for i in range(1,4):
				score_R = aligner.score(memory_pattern_list_R[i], potential_R)
				memory_score_list.append(score_R)
			del memory_score_list[1:4]
			distance_score = float(sum(memory_distance_list))/len(memory_distance_list)
			alignment_score = float(sum(memory_score_list))/len(memory_score_list)
			alignment_score1 = alignment_score/len(potential_R)
			repeat_score = (float(sum(memory_score_list)-sum(memory_distance_list))/(repeat_end-repeat_start+1))*repeat_time
			if repeat_score >= alignAS*repeat_time and distance_score <= gapAS: 
				return potential_R.upper(),repeat_time,repeat_score,len(potential_R.upper()),repeat_start+based,repeat_end+based,variant_start,alignment_score1,distance_score


def count_microsatellite(pattern_start, potential_R, search_size, ref, variant_start,based):
	repeat_time = 1
	repeat_start = pattern_start
	repeat_end = pattern_start + search_size-1

	repeat_start1 = repeat_start
	repeat_end1 = repeat_end

	pattern_start1 = pattern_start

	potential_R = Seq(potential_R.upper())
	while potential_R == Seq(ref[pattern_start+search_size:pattern_start+search_size+search_size].upper()):
		pattern_start = pattern_start+search_size
		repeat_time = repeat_time +1
		repeat_end = repeat_end + search_size
	while potential_R == Seq(ref[pattern_start1-search_size:pattern_start1].upper()):
		pattern_start1 = pattern_start1 - search_size
		repeat_time = repeat_time +1
		repeat_start = repeat_start -search_size
	if repeat_start != repeat_start1 or repeat_end != repeat_end1:
		if repeat_time >=5:
			distance_score = 0
			alignment_score = 5
			repeat_score = repeat_time*5 
			return potential_R.upper(),repeat_time,repeat_score,len(potential_R.upper()),repeat_start+based,repeat_end+based,variant_start,alignment_score,distance_score












































