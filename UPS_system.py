import read_reference

def left_coordination(start_C,pattern_C,pattern_D): 
	if len(pattern_C)>len(pattern_D):
		N_start_C=int(start_C)+1 
	elif len(pattern_C)<len(pattern_D):
		N_start_C=int(start_C)   

	return N_start_C

def right_coordination(start_C,pattern_C,pattern_D,b,based): 
	if len(pattern_C)>len(pattern_D):
		variants_C = pattern_C.upper()[1:] 
		N_end_C = int(start_C)+len(pattern_C)-based
		flank_right = b[N_end_C].upper() 
	elif len(pattern_C)<len(pattern_D):
		variants_C = pattern_D.upper()[1:]
		N_end_C = int(start_C)-based+1
		flank_right = b[N_end_C].upper()

	while variants_C[0]==flank_right:
		N_end_C=N_end_C+1
		flank_right = b[N_end_C].upper()
		variants_C=variants_C[1:]+variants_C[0]

	N_end_C=N_end_C-1+based

	return N_end_C  

def UPS(df,ref_file,based):
	LEFT=[]
	RIGHT=[]
	RIGHT_align=[]
	LEFT_2=[]
	RIGHT_2=[]
	chr_name=""
	for i in range(0,df.shape[0]):
		chromosome=df.loc[i,"Chromosome"]
		if chr_name!=chromosome:
			chr_name=chromosome
			ref=read_reference.read_reference_sequence(ref_file,chr_name)
		
		if df.loc[i,"REF"].upper()[0]==df.loc[i,"ALT"].upper()[0]:
			variants1 = left_coordination(df.loc[i,"Position"],df.loc[i,"REF"],df.loc[i,"ALT"])
			LEFT.append(variants1)
			variants2 = right_coordination(df.loc[i,"Position"],df.loc[i,"REF"],df.loc[i,"ALT"],ref,based)
			RIGHT.append(variants2)

			if len(df.loc[i,"REF"])>len(df.loc[i,"ALT"]):
				R_P=variants2-(len(df.loc[i,"REF"])-len(df.loc[i,"ALT"]))+1
				RIGHT_align.append(R_P)
				if R_P>(df.loc[i,"Position"]+1):
					LEFT_2.append(variants1)
					RIGHT_2.append(variants2)
				else:
					LEFT_2.append("")
					RIGHT_2.append("")
				
			elif len(df.loc[i,"REF"])<len(df.loc[i,"ALT"]):
				RIGHT_align.append(variants2)
				if variants2>df.loc[i,"Position"]:
					LEFT_2.append(variants1)
					RIGHT_2.append(variants2)
				else:
					LEFT_2.append("")
					RIGHT_2.append("")
			
		elif df.loc[i,"REF"].upper()[0]!=df.loc[i,"ALT"].upper()[0]:
			LEFT.append(int(df.loc[i,"Position"]))
			RIGHT.append(int(df.loc[i,"Position"])+len(df.loc[i,"REF"])-1)
			RIGHT_align.append(int(df.loc[i,"Position"]))
			LEFT_2.append("")
			RIGHT_2.append("")

	return LEFT,RIGHT,RIGHT_align,LEFT_2,RIGHT_2





















