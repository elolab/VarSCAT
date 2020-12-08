from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import read_reference

def mutation_seq(df,ref_file):
    REF_seq=[]
    MUT_seq=[]
    chr_name=""
    for i in range(0,df.shape[0]):
        # read reference
        chromosome=df.loc[i,"Chromosome"]
        if chr_name!=chromosome:
            chr_name=chromosome
            ref=read_reference.read_reference_sequence(ref_file,chr_name)
        # Insertion UPS_LEFT is reference position, while deletion UPS_LEFT is the first variants position
        if len(df.loc[i,"REF"]) > len(df.loc[i,"ALT"]):
            # -1-1, means first -1: convert Reference genome 1-based to 0-based, second -1: output reference position 
            # for end position, 1-based to 0-based naturalized python last coordinate not included
                
        # "UPS_LEFT": first variants position in 1-based, -1, to reference position, -1 to one base left to variants position, -1 convert to 0-based
            S = df.loc[i,"UPS_L"]-3
        # "UPS_RIGHT": last variants position in 1-based, +1, to one base right to last variants position, -1 convert to 0-based, 
            E = df.loc[i,"UPS_R"]+1
            Ref_seq = Seq(ref[S:E].upper())
            Mut_seq = Seq(ref[S].upper())+Seq(df.loc[i,"ALT"][0:].upper())+" -"*(len(df.loc[i,"REF"])-1)+Seq(ref[df.loc[i,"UPS_L"]+len(df.loc[i,"REF"])-2:E].upper())
            REF_seq.append(str(Ref_seq))
            MUT_seq.append(str(Mut_seq))
        elif len(df.loc[i,"REF"]) < len(df.loc[i,"ALT"]):
        # "UPS_LEFT": reference position in 1-based,-1 to one base left to variants position, -1 convert to 0-based
            S = df.loc[i,"UPS_L"]-2
        # "UPS_RIGHT": last UPS position in reference in 1-based, +1, to one base right to last variants position, -1 convert to 0-based, +1, python last corrdinate not included
            E = df.loc[i,"UPS_R"]+1
            Ref_seq = Seq(ref[S:E].upper())
            Mut_seq = Seq(ref[S].upper())+Seq(df.loc[i,"ALT"].upper())+Seq(ref[df.loc[i,"UPS_L"]:E].upper())                
            REF_seq.append(str(Ref_seq))
            MUT_seq.append(str(Mut_seq))
        elif len(df.loc[i,"REF"]) == len(df.loc[i,"ALT"]): # no multi-alts marked with ",", due to it converted to single alt in read_VCF.py
            S = df.loc[i,"UPS_L"]-2
            E = df.loc[i,"UPS_R"]+1
            Ref_seq = Seq(ref[S:E].upper())
            Mut_seq = Seq(ref[S].upper())+df.loc[i,"ALT"]+Seq(ref[df.loc[i,"UPS_L"]+len(df.loc[i,"REF"])-1:E].upper())
            REF_seq.append(str(Ref_seq))
            MUT_seq.append(str(Mut_seq))
    
    return REF_seq, MUT_seq






















































