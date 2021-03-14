from ordered_set import OrderedSet

def combine_complex(df):
	chr_set=OrderedSet(df.loc[:,"Chromosome"])
	CLS_all=[]
	distance_to_3T_all=[]
	for chr_name in chr_set:
		df_chr=df.loc[df["Chromosome"]==chr_name].copy()
		df_chr.reset_index(drop=True,inplace=True)	
		
		left_coordinates = []
		right_coordinates = []
		for i in range(0,df_chr.shape[0]):
			left_coordinates.append(int(df_chr.loc[i,"UPS_L"]))
			right_coordinates.append(int(df_chr.loc[i,"UPS_R"]))

		left_coordinates = left_coordinates[1:]
		right_coordinates = right_coordinates[:len(right_coordinates)-1]
	
		k=0
		m=0
		n=0
		e=0
		f=0
		CLS=[]
		distance_to_3T = []
		for i in range(0,len(left_coordinates)):
			distance=left_coordinates[i]-right_coordinates[i]
			distance_to_3T.append(distance)

			if df_chr.loc[i+1,"REF"]>=df_chr.loc[i+1,"ALT"]:
				criteria=1
			elif df_chr.loc[i+1,"REF"]<df_chr.loc[i+1,"ALT"]:
				criteria=0

			if distance<=criteria: 
				if m!=n:
					k=k+1
				CLS1 = df_chr.loc[i,"Chromosome"]+"_"+df_chr.loc[i,"Position"]+"_"+df_chr.loc[i,"REF"]+"_"+df_chr.loc[i,"ALT"]+"_"+df_chr.loc[i,"Genotype"]+"_"+"ADJ"+"_"+str(k)
				CLS2 = df_chr.loc[i+1,"Chromosome"]+"_"+df_chr.loc[i+1,"Position"]+"_"+df_chr.loc[i+1,"REF"]+"_"+df_chr.loc[i+1,"ALT"]+"_"+df_chr.loc[i+1,"Genotype"]+"_"+"ADJ"+"_"+str(k)
				if CLS1 not in CLS:
					CLS.append(CLS1)
				if CLS2 not in CLS:    
					CLS.append(CLS2)
					f=f+1
				m=n
			else:
				if e!=f:
					e=f
					n=n+1
				else:
					CLS.append("")
					e=f
					n=n+1

		CLS.append("")
		distance_to_3T.append("")
		
		if len(CLS)!=len(distance_to_3T):
			del CLS[-1]
	
		CLS_all=CLS_all+CLS
		distance_to_3T_all=distance_to_3T_all+distance_to_3T	

	return CLS_all,distance_to_3T_all
	






























































