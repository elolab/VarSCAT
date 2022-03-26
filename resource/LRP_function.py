from .read_reference import read_reference_sequence

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

def LRP(df,ref_file,based):
	LEFT=[]
	RIGHT=[]
	RIGHT_align=[]
	chr_name=""
	for i in range(0,len(df)):
		chromosome=df[i][0]
		if chr_name!=chromosome:
			chr_name=chromosome
			ref=read_reference_sequence(ref_file,chr_name)
		
		if df[i][2].upper()[0]==df[i][3].upper()[0]:
			variants1 = left_coordination(df[i][1],df[i][2],df[i][3])
			LEFT.append(variants1)
			variants2 = right_coordination(df[i][1],df[i][2],df[i][3],ref,based)
			RIGHT.append(variants2)

			if len(df[i][2])>len(df[i][3]):
				R_P=variants2-(len(df[i][2])-len(df[i][3]))+1
				RIGHT_align.append(R_P)
				
			elif len(df[i][2])<len(df[i][3]):
				RIGHT_align.append(variants2)
			
		elif df[i][2].upper()[0]!=df[i][3].upper()[0]:
			LEFT.append(int(df[i][1]))
			RIGHT.append(int(df[i][1])+len(df[i][2])-1)
			RIGHT_align.append(int(df[i][1]))

	return LEFT,RIGHT_align,RIGHT
