from Bio.Seq import Seq
import read_reference

def mutation_seq(df,ref_file,based):
	REF_seq=[]
	MUT_seq=[]
	chr_name=""
	for i in range(0,df.shape[0]):
		chromosome=df.loc[i,"Chromosome"]
		if chr_name!=chromosome:
			chr_name=chromosome
			ref=read_reference.read_reference_sequence(ref_file,chr_name)
		if len(df.loc[i,"REF"]) > len(df.loc[i,"ALT"]):
			
			if df.loc[i,"REF"].upper()[0]==df.loc[i,"ALT"].upper()[0]:
				S = df.loc[i,"UPS_L"]-1-based
				E = df.loc[i,"UPS_R"]+1+1-based
				Ref_seq = Seq(ref[S:E].upper())
				Mut_seq = Seq(ref[S].upper())+" -"*(len(df.loc[i,"REF"])-1)+Seq(ref[df.loc[i,"UPS_L"]+len(df.loc[i,"REF"])-1-based:E].upper())			
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
			elif df.loc[i,"REF"].upper()[0]!=df.loc[i,"ALT"].upper()[0]:
				S = df.loc[i,"UPS_L"]-1-based
				E = df.loc[i,"UPS_R"]+1+1-based
				Ref_seq = Seq(ref[S:E].upper())
				Mut_seq = Seq(ref[S].upper())+Seq(df.loc[i,"ALT"][0:].upper())+" -"*(len(df.loc[i,"REF"])-1)+Seq(ref[df.loc[i,"UPS_R"]+len(df.loc[i,"REF"])-1-based:E].upper())			
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
		elif len(df.loc[i,"REF"]) < len(df.loc[i,"ALT"]):
			if df.loc[i,"REF"].upper()[0]==df.loc[i,"ALT"].upper()[0]:
				S = df.loc[i,"UPS_L"]-based 
				E = df.loc[i,"UPS_R"]+1+1-based 
				Ref_seq = Seq(ref[S:E].upper())

				Mut_seq = Seq(df.loc[i,"ALT"].upper())+Seq(ref[S+1:E].upper())				
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
			elif df.loc[i,"REF"].upper()[0]!=df.loc[i,"ALT"].upper()[0]:
				S = df.loc[i,"UPS_L"]-1-based 
				E = df.loc[i,"UPS_R"]+1+1-based 
				Ref_seq = Seq(ref[S:E].upper())
				Mut_seq = Seq(ref[S].upper())+Seq(df.loc[i,"ALT"].upper())+Seq(ref[S+1+1:E].upper())
				REF_seq.append(str(Ref_seq))
				MUT_seq.append(str(Mut_seq))
		elif len(df.loc[i,"REF"]) == len(df.loc[i,"ALT"]): 
			S = df.loc[i,"UPS_L"]-1-based
			E = df.loc[i,"UPS_R"]+1+1-based 
			Ref_seq = Seq(ref[S:E].upper())
			Mut_seq = Seq(ref[S].upper())+df.loc[i,"ALT"]+Seq(ref[df.loc[i,"UPS_L"]+len(df.loc[i,"REF"])-based:E].upper())
			REF_seq.append(str(Ref_seq))
			MUT_seq.append(str(Mut_seq))
	
	return REF_seq, MUT_seq






















































