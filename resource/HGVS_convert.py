from Bio.Seq import Seq
from .read_reference import read_reference_sequence
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
	for i in range(0,len(df)):
		chromosome=df[i][0]
		if chr_name!=chromosome:
			chr_name=chromosome
			ref=read_reference_sequence(ref_file,chr_name)

		
		g_prefix="g"
		if chromosome=="MT" or chromosome=="chrMT":
			g_prefix="m"

		if df[i][2].upper()[0]==df[i][3].upper()[0]:
			if len(df[i][2])>len(df[i][3]):
				if int(df[i][5])==int(df[i][1])+1:
					if len(df[i][2])-len(df[i][3])==1:
						item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"del"
					elif len(df[i][2])-len(df[i][3])>1:
						item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"_"+str(int(df[i][5])+len(df[i][2])-len(df[i][3])-1)+"del"
				elif int(df[i][5])!=int(df[i][1])+1:
					if len(df[i][2])-len(df[i][3])>1:
						Right_align_seq=ref[int(df[i][5])-based:int(df[i][6])-based+1].upper()
						in_seq_result=in_seq_pattern(Right_align_seq)
						in_seq_RT=in_seq_result[0]
						in_seq_RM=in_seq_result[1]
						
						out_seq_result=out_seq_pattern(in_seq_RM,ref,based,int(df[i][5]))
						out_seq_left=out_seq_result[0]
						out_seq_right=out_seq_result[1]
						out_seq_RT=out_seq_result[2]
		
						if out_seq_RT==0:
							item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"_"+str(int(df[i][5])+len(df[i][2])-len(df[i][3])-1)+"del"
						elif out_seq_RT!=0:
							if len(in_seq_RM)==1:
								item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"_"+str(int(df[i][5])+len(df[i][2])-len(df[i][3])-1)+"del"
							elif len(in_seq_RM)>1:
								item = str(df[i][0])+":"+g_prefix+"."+str(out_seq_left)+in_seq_RM+"["+str(out_seq_RT)+"]"						
					elif len(df[i][2])-len(df[i][3])==1:
						item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"del"

			elif len(df[i][2])<len(df[i][3]):
				if int(df[i][5])==int(df[i][1]):
					item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"_"+str(df[i][5]+1)+"ins"+str(df[i][3][1:])
				elif int(df[i][5])!=int(df[i][1]):
					if len(df[i][3])-len(df[i][2])==1:
						item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"dup"
					elif len(df[i][3])-len(df[i][2])>1:
						Right_align_seq=form_shifted_variant(df[i][3][1:],int(df[i][5]),int(df[i][1]))
						in_seq_result=in_seq_pattern(Right_align_seq)
						in_seq_RT=in_seq_result[0]
						in_seq_RM=in_seq_result[1]
						
						out_seq_result=out_seq_pattern(in_seq_RM,ref,based,int(df[i][5])+1) 
						out_seq_left=out_seq_result[0]
						out_seq_right=out_seq_result[1]
						out_seq_RT=out_seq_result[2]
						if out_seq_RT==0:
							item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"_"+str(df[i][5]+1)+"ins"+Right_align_seq
						elif out_seq_RT!=0:
							if len(in_seq_RM)==1:
								if out_seq_RT<len(df[i][3][1:]):
									item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5])+"_"+str(df[i][5]+1)+"ins"+str(df[i][3][1:])
								elif out_seq_RT>=len(df[i][3][1:]):
									item = str(df[i][0])+":"+g_prefix+"."+str(df[i][5]-len(df[i][3])+1+1)+"_"+str(df[i][5])+"dup"
							elif len(in_seq_RM)>1:
								if (out_seq_RT+in_seq_RT)==2: 
									item = str(df[i][0])+":"+g_prefix+"."+str(out_seq_left)+"_"+str(out_seq_right)+"dup"
								elif (out_seq_RT+in_seq_RT)>2:
									item = str(df[i][0])+":"+g_prefix+"."+str(out_seq_left)+in_seq_RM+"["+str(out_seq_RT+in_seq_RT)+"]"		

		elif df[i][2].upper()[0]!=df[i][3].upper()[0]:
			if len(df[i][2])==1 and len(df[i][3])==1:
				item = str(df[i][0])+":"+g_prefix+"."+str(df[i][1])+df[i][2]+">"+df[i][3]
			elif len(df[i][2])!=1 or len(df[i][3])!=1:	 			
				if len(df[i][2])==1:
					item = str(df[i][0])+":"+g_prefix+"."+str(df[i][1])+"delins"+df[i][3]
				elif len(df[i][2])!=1 and len(df[i][2])!=len(df[i][3]):
					item = str(df[i][0])+":"+g_prefix+"."+str(df[i][1])+"_"+str(int(df[i][1])+len(df[i][2])-1)+"delins"+df[i][3]
				elif len(df[i][2])!=1 and len(df[i][2])==len(df[i][3]):
					if Seq(df[i][2]).reverse_complement()==df[i][3]:
						item = str(df[i][0])+":"+g_prefix+"."+str(df[i][1])+"_"+str(int(df[i][1])+len(df[i][2])-1)+"inv"
					else:
						item = str(df[i][0])+":"+g_prefix+"."+str(df[i][1])+"_"+str(int(df[i][1])+len(df[i][2])-1)+"delins"+df[i][3]

		HGVS_result.append(item)

	return HGVS_result

