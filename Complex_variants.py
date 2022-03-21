from ordered_set import OrderedSet

def combine_complex(df):
	chr_set=OrderedSet(df.loc[:,"Chromosome"])
	distance_to_3T_all=[]
	for chr_name in chr_set:
		df_chr=df.loc[df["Chromosome"]==chr_name].copy()
		df_chr.reset_index(drop=True,inplace=True)	
		
		left_coordinates = []
		right_coordinates = []
		for i in range(0,df_chr.shape[0]):
			left_coordinates.append(int(df_chr.loc[i,"5'_aligned"]))
			right_coordinates.append(int(df_chr.loc[i,"3'_edge"]))

		left_coordinates = left_coordinates[1:]
		right_coordinates = right_coordinates[:len(right_coordinates)-1]
	
		distance_to_3T = []
		for i in range(0,len(left_coordinates)):
			distance=left_coordinates[i]-right_coordinates[i]
			if distance<0:
				distance=0
			distance_to_3T.append(distance)

		distance_to_3T.append("")
		
		distance_to_3T_all=distance_to_3T_all+distance_to_3T

	return distance_to_3T_all

