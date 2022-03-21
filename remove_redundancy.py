from Bio.Seq import Seq
from Bio import Align
from Bio import motifs
from Bio.SeqUtils import GC

def ATGC_diff(seq,pattern):
	seq_A=seq.count("A")/len(seq)
	seq_T=seq.count("T")/len(seq)
	seq_G=seq.count("G")/len(seq)
	seq_C=seq.count("C")/len(seq)
	pattern_A=pattern.count("A")/len(pattern)
	pattern_T=pattern.count("T")/len(pattern)
	pattern_G=pattern.count("G")/len(pattern)
	pattern_C=pattern.count("C")/len(pattern)

	diff=abs(seq_A-pattern_A)+abs(seq_T-pattern_T)+abs(seq_G-pattern_G)+abs(seq_C-pattern_C)
	pattern_freq=[pattern_A,pattern_T,pattern_G,pattern_C]	
	
	return diff,pattern_freq
	
def clean_df(df,ref,base): 
	df = df.sort_values('RepeatS',ascending=False) 
	df = df.reset_index()
	drop_list = []
	for i in range(0,df.shape[0]-1): 
		for j in range(i+1,df.shape[0]):
			region_ij_min = float(max(df.iloc[i]['End'],df.iloc[j]['End'])-min(df.iloc[i]['Start'],df.iloc[j]['Start'])+1)
			region_ij_max = float(min(df.iloc[i]['End'],df.iloc[j]['End'])-max(df.iloc[i]['Start'],df.iloc[j]['Start'])+1)
			redundant_score = float(region_ij_max/region_ij_min)
			if redundant_score > 0.5:
				diff_i=ATGC_diff(ref[int(df.iloc[i]['Start']-base):int(df.iloc[i]['End']-base+1)].upper(),df.iloc[i]['Motifs'])
				diff_j=ATGC_diff(ref[int(df.iloc[j]['Start']-base):int(df.iloc[j]['End']-base+1)].upper(),df.iloc[j]['Motifs'])
				if diff_i[1]==diff_j[1]:
					drop_list.append(j)
				elif diff_i[1]!=diff_j[1] and diff_i[0]<=diff_j[0]:
					drop_list.append(j)
				elif diff_i[1]!=diff_j[1] and diff_i[0]>diff_j[0]:
					drop_list.append(i)
					break

	drop_list = set(drop_list)
	drop_list = list(drop_list)
	df = df.drop(df.index[drop_list])

	return df

