from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import read_reference
import math

def form_shifted_variant(seq,RA,OP):
	epoch=RA-OP
	while epoch>0:
		seq=seq[1:]+seq[0]
		epoch=epoch-1
	
	return seq
		
def FactorCalu(num):
	divisor_num=[]
	count=num//2 
	while count>0:
		if num%count==0:
			divisor_num.append(count)
			count=count-1
		else:	
			count=count-1

	divisor_num.append(num)

	return divisor_num

def in_seq_pattern(variant_seq):
	divisor=FactorCalu(len(variant_seq))
	repeat_time=[]
	repeat_motif=[]
	for m in divisor:
		b=len(variant_seq)
		a=len(variant_seq)-m
		pattern=variant_seq[a:b]
		RT=1 
		RM=pattern 
		while a>0:
			a=a-m
			b=b-m
			pattern_next=variant_seq[a:b]
			if pattern==pattern_next:
				RT=RT+1
			else:
				break
				
		if RT*m==len(variant_seq):
			repeat_time.append(RT)
			repeat_motif.append(pattern)
	
	if len(repeat_time)>1:
		repeat_time2=repeat_time[-2]
		repeat_motif2=repeat_motif[-2]
	else:
		repeat_time2=repeat_time[-1]
		repeat_motif2=repeat_motif[-1]
			

	return repeat_time2,repeat_motif2
			

def out_seq_pattern(pattern,reference,based,position):
	a=position-based-len(pattern)
	b=position-based
	pattern_next=reference[a:b].upper()
	RT=0 
	while pattern==pattern_next:
		RT=RT+1
		a=a-len(pattern)
		b=b-len(pattern)
		pattern_next=reference[a:b].upper()

	a=a+len(pattern)+based
	b=b+len(pattern)-1+based 

	return a,b,RT

def HGVS_transform(df,ref_file,based):
	HGVS_result=[]
	chr_name=""
	for i in range(0,df.shape[0]):
		chromosome=df.loc[i,"Chromosome"]
		if chr_name!=chromosome:
			chr_name=chromosome
			ref=read_reference.read_reference_sequence(ref_file,chr_name)
		if df.loc[i,"REF"].upper()[0]==df.loc[i,"ALT"].upper()[0]:
			if len(df.loc[i,"REF"])>len(df.loc[i,"ALT"]):
				if int(df.loc[i,"Right align position"])==int(df.loc[i,"Position"])+1:
					if len(df.loc[i,"REF"])-len(df.loc[i,"ALT"])==1:
						item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"del"
					elif len(df.loc[i,"REF"])-len(df.loc[i,"ALT"])>1:
						item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"_"+str(int(df.loc[i,"Right align position"])+len(df.loc[i,"REF"])-len(df.loc[i,"ALT"])-1)+"del"
				elif int(df.loc[i,"Right align position"])!=int(df.loc[i,"Position"])+1:
					if len(df.loc[i,"REF"])-len(df.loc[i,"ALT"])>1:
						Right_align_seq=ref[int(df.loc[i,"Right align position"])-based:int(df.loc[i,"UPS_R"])-based+1].upper()
						in_seq_result=in_seq_pattern(Right_align_seq)
						in_seq_RT=in_seq_result[0]
						in_seq_RM=in_seq_result[1]
						
						out_seq_result=out_seq_pattern(in_seq_RM,ref,based,int(df.loc[i,"Right align position"]))
						out_seq_left=out_seq_result[0]
						out_seq_right=out_seq_result[1]
						out_seq_RT=out_seq_result[2]
		
						if out_seq_RT==0:
							item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"_"+str(int(df.loc[i,"Right align position"])+len(df.loc[i,"REF"])-len(df.loc[i,"ALT"])-1)+"del"
						elif out_seq_RT!=0:
							if len(in_seq_RM)==1:
								item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"_"+str(int(df.loc[i,"Right align position"])+len(df.loc[i,"REF"])-len(df.loc[i,"ALT"])-1)+"del"
							elif len(in_seq_RM)>1:
								item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(out_seq_left)+in_seq_RM+"["+str(out_seq_RT)+"]"						
					elif len(df.loc[i,"REF"])-len(df.loc[i,"ALT"])==1:
						item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"del"

			elif len(df.loc[i,"REF"])<len(df.loc[i,"ALT"]):
				if int(df.loc[i,"Right align position"])==int(df.loc[i,"Position"]):
					item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"_"+str(df.loc[i,"Right align position"]+1)+"ins"+str(df.loc[i,"ALT"][1:])
				elif int(df.loc[i,"Right align position"])!=int(df.loc[i,"Position"]):
					if len(df.loc[i,"ALT"])-len(df.loc[i,"REF"])==1:
						item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"dup"
					elif len(df.loc[i,"ALT"])-len(df.loc[i,"REF"])>1:
						Right_align_seq=form_shifted_variant(df.loc[i,"ALT"][1:],int(df.loc[i,"Right align position"]),int(df.loc[i,"Position"]))
						in_seq_result=in_seq_pattern(Right_align_seq)
						in_seq_RT=in_seq_result[0]
						in_seq_RM=in_seq_result[1]
						
						out_seq_result=out_seq_pattern(in_seq_RM,ref,based,int(df.loc[i,"Right align position"])+1) 
						out_seq_left=out_seq_result[0]
						out_seq_right=out_seq_result[1]
						out_seq_RT=out_seq_result[2]
						if out_seq_RT==0:
							item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"_"+str(df.loc[i,"Right align position"]+1)+"ins"+in_seq_RM
						elif out_seq_RT!=0:
							if len(in_seq_RM)==1:
								if out_seq_RT<len(df.loc[i,"ALT"][1:]):
									item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"])+"_"+str(df.loc[i,"Right align position"]+1)+"ins"+str(df.loc[i,"ALT"][1:])
								elif out_seq_RT>=len(df.loc[i,"ALT"][1:]):
									item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Right align position"]-len(df.loc[i,"ALT"])+1+1)+"_"+str(df.loc[i,"Right align position"])+"dup"
							elif len(in_seq_RM)>1:
								if out_seq_RT==1: 
									item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(out_seq_left)+"_"+str(out_seq_right)+"dup"
								elif out_seq_RT>1:
									item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(out_seq_left)+in_seq_RM+"["+str(out_seq_RT+in_seq_RT)+"]"		

		elif df.loc[i,"REF"].upper()[0]!=df.loc[i,"ALT"].upper()[0]:
			if len(df.loc[i,"REF"])==1 and len(df.loc[i,"ALT"])==1:
				item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Position"])+df.loc[i,"REF"]+">"+df.loc[i,"ALT"] 
			elif len(df.loc[i,"REF"])!=1 or len(df.loc[i,"ALT"])!=1:	 			
				if len(df.loc[i,"REF"])==1:
					item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Position"])+"delins"+df.loc[i,"ALT"]
				elif len(df.loc[i,"REF"])!=1 and len(df.loc[i,"REF"])!=len(df.loc[i,"ALT"]):
					item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Position"])+"_"+str(int(df.loc[i,"Position"])+len(df.loc[i,"REF"])-1)+"delins"+df.loc[i,"ALT"] 
				elif len(df.loc[i,"REF"])!=1 and len(df.loc[i,"REF"])==len(df.loc[i,"ALT"]):
					if Seq(df.loc[i,"REF"]).reverse_complement()==df.loc[i,"ALT"]:
						item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Position"])+"_"+str(int(df.loc[i,"Position"])+len(df.loc[i,"REF"])-1)+"inv"
					else:
						item = str(df.loc[i,"Chromosome"])+":"+"g"+"."+str(df.loc[i,"Position"])+"_"+str(int(df.loc[i,"Position"])+len(df.loc[i,"REF"])-1)+"delins"+df.loc[i,"ALT"]

		HGVS_result.append(item)

	return HGVS_result










						

