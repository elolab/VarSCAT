from .read_reference import read_reference_sequence

def mutation_seq(df,ref_file,based):
	REF_seq=[]
	MUT_seq=[]
	chr_name=""
	for i in range(0,len(df)):
		chromosome=df[i][0]
		if chr_name!=chromosome:
			chr_name=chromosome
			ref=read_reference_sequence(ref_file,chr_name)
		if len(df[i][2]) > len(df[i][3]):
			
			if df[i][2].upper()[0]==df[i][3].upper()[0]:
				S = df[i][4]-1-based
				E = df[i][6]+1+1-based
				Ref_seq = ref[S:E].upper()
				Mut_seq = ref[S].upper()+"-"*(len(df[i][2])-1)+ref[df[i][4]+len(df[i][2])-1-based:E].upper()	
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
			elif df[i][2].upper()[0]!=df[i][3].upper()[0]:
				S = df[i][4]-1-based
				E = df[i][6]+1+1-based
				Ref_seq = ref[S:E].upper()
				Mut_seq = ref[S].upper()+df[i][3][0:].upper()+"-"*(len(df[i][2])-1)+ref[df[i][6]+len(df[i][2])-1-based:E].upper()			
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
		elif len(df[i][2]) < len(df[i][3]):
			if df[i][2].upper()[0]==df[i][3].upper()[0]:
				S = df[i][4]-based 
				E = df[i][6]+1+1-based 
				Ref_seq = ref[S:E].upper()

				Mut_seq = df[i][3].upper()+ref[S+1:E].upper()		
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
			elif df[i][2].upper()[0]!=df[i][3].upper()[0]:
				S = df[i][4]-1-based 
				E = df[i][6]+1+1-based 
				Ref_seq = ref[S:E].upper()
				Mut_seq = ref[S].upper()+df[i][3].upper()+ref[S+1+1:E].upper()
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
		elif len(df[i][2]) == len(df[i][3]): 
			S = df[i][4]-1-based
			E = df[i][6]+1+1-based 
			Ref_seq = ref[S:E].upper()
			Mut_seq = ref[S].upper()+df[i][3]+ref[df[i][4]+len(df[i][2])-based:E].upper()
			REF_seq.append(str(Ref_seq))
			MUT_seq.append(str(Mut_seq))
	
	return REF_seq, MUT_seq
