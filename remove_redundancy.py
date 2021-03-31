from Bio.Seq import Seq
from Bio import Align
from Bio import motifs

def clean_df(df): 
	df = df.sort_values('RepeatS',ascending=False) 
	df = df.reset_index()
	drop_list = []
	for i in range(0,df.shape[0]-1): 
		for j in range(i+1,df.shape[0]):
			region_ij = float(max(df.iloc[i]['End'],df.iloc[j]['End'])-min(df.iloc[i]['Start'],df.iloc[j]['Start'])+1)
			redundant_score = float(min(df.iloc[i]['End'],df.iloc[j]['End'])-max(df.iloc[i]['Start'],df.iloc[j]['Start'])+1)/region_ij
			if redundant_score > 0.3:
				drop_list.append(j)

	drop_list = set(drop_list)
	drop_list = list(drop_list)
	df = df.drop(df.index[drop_list])

	return df
