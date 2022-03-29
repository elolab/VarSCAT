from ordered_set import OrderedSet

def combine_complex(df):
	chr_set=OrderedSet(row[0] for row in df)
	distance_to_3T_all=[]
	for chr_name in chr_set:
		df_chr = [x for x in df if x[0] == chr_name].copy()	
		
		left_coordinates = []
		right_coordinates = []
		for i in range(0,len(df_chr)):
			left_coordinates.append(int(df_chr[i][4]))
			right_coordinates.append(int(df_chr[i][6]))

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
