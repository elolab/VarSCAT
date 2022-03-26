from operator import itemgetter

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
	df = sorted(df, key=itemgetter(5),reverse=True)
	drop_list = []
	for i in range(0,len(df)-1): 
		for j in range(i+1,len(df)):
			region_ij_min = float(max(df[i][4],df[j][4])-min(df[i][3],df[j][3])+1)
			region_ij_max = float(min(df[i][4],df[j][4])-max(df[i][3],df[j][3])+1)
			redundant_score = float(region_ij_max/region_ij_min)
			if redundant_score > 0.5:
				diff_i=ATGC_diff(ref[int(df[i][3]-base):int(df[i][4]-base+1)].upper(),df[i][0])
				diff_j=ATGC_diff(ref[int(df[j][3]-base):int(df[j][4]-base+1)].upper(),df[j][0])
				if diff_i[1]==diff_j[1]:
					drop_list.append(j)
				elif diff_i[1]!=diff_j[1] and diff_i[0]<=diff_j[0]:
					drop_list.append(j)
				elif diff_i[1]!=diff_j[1] and diff_i[0]>diff_j[0]:
					drop_list.append(i)
					break

	drop_list = set(drop_list)
	drop_list = list(drop_list)
	for index in sorted(drop_list, reverse=True):
		del df[index]

	return df
