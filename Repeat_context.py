import sys
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Align
from Bio import motifs
import repeatAlign # strict tandem repeat (microsatellite function, for 1-5bp) & mismatch-tolerate tandem repeat (minisatellite function)
import remove_redundancy
import Local_GC
import read_reference

# kick out the non-significant repeat, redundancy
#def select_significant(pattern_start, potential_R, search_size,ref):


# loop 1: catching variant
# VCF 1-based, python 0-based, So REF start at position-1, but variant is ref+variant (TGG  G), variant position should be position+1, turn to be position-1+1
# this variant_start is python 0-based
def catching_variant_Indel(position_l1, REF_l1, ALT_l1, Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1):
    variant_start = int(position_l1)-1 # 18.11: due to there might be complex variant (TGG  A), so the reference position need to be considered
    if len(REF_l1) > len(ALT_l1):
        variant_pattern = ALT_l1[0]+REF_l1[1:]
        variant_length = len(REF_l1)
    elif len(REF_l1) < len(ALT_l1):
        variant_pattern = ALT_l1
        variant_length = len(ALT_l1)
    Repeat_raw_L1 = sizing_window(variant_start,variant_pattern,variant_length,Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1) # 1 here is the min length of repeat
    return Repeat_raw_L1

def catching_variant_SNV(position_l1, REF_l1, ALT_l1, Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1):
    variant_start = int(position_l1)-1
    variant_pattern = REF_l1
    variant_length = len(REF_l1)
    Repeat_raw_L1 = sizing_window(variant_start,variant_pattern,variant_length,Ref1,min_size1,max_size1,exact_size1,alignAS1,gapAS1) # 1 here is the min length of repeat
    return Repeat_raw_L1

# loop 2: choose different sizes of window i
def sizing_window(start,pattern,length, Ref2,min_size2,max_size2,exact_size2,alignAS2,gapAS2):
    Repeat_df_L2 = pd.DataFrame(columns=["Motifs", "Period","RepeatS","Size","Start","End","Position_O","AlignmentS","GapS"])
    ###############################################
    for i in range(min_size2,max_size2): # I set max_size=10, min_size=1 (marked as 11, due to python not included last) as default, 
        #sliding_window_micro(i,start,length)
        if i <=exact_size2: # here is the max length of microsatellite, user-defined
            Repeat_raw_L2 = sliding_window_micro(i,start,length,Ref2)
            if Repeat_raw_L2.shape[0]>1:
                Repeat_raw_L2_sort = remove_redundancy.clean_df(Repeat_raw_L2)
                Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2_sort,ignore_index = True)
            elif Repeat_raw_L2.shape[0]==1:
                Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2,ignore_index = True)
        elif i >exact_size2:
            Repeat_raw_L2 = sliding_window_mini(i,start,length,Ref2,alignAS2,gapAS2)
            if Repeat_raw_L2.shape[0]>1:
                Repeat_raw_L2_sort = remove_redundancy.clean_df(Repeat_raw_L2)
                Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2_sort,ignore_index = True)
            elif Repeat_raw_L2.shape[0]==1:
                Repeat_df_L2=Repeat_df_L2.append(Repeat_raw_L2,ignore_index = True)

    # at this stage should kick redundancy
    if Repeat_df_L2.shape[0]>1:
        Repeat_df_L2=remove_redundancy.clean_df(Repeat_df_L2)
    return Repeat_df_L2

# loop 3: sliding the certain size of window to different positions j
def sliding_window_micro(window_size,V_start,V_length,b):
    # all the possible windows at least one base overlap with variants
    Repeat_df_L3 = pd.DataFrame(columns=["Motifs", "Period","RepeatS","Size","Start","End","Position_O","AlignmentS","GapS"])
    for j in range(V_start-window_size+1,V_start+V_length): # the last window position should be V_start+V_length-1, but python do not include last index, so V_start+V_length
        V_pattern = b[j:j+window_size]
        Repeat_raw_L3 = repeatAlign.count_microsatellite(j,V_pattern,window_size,b,V_start)
        if Repeat_raw_L3 != None:
            Repeat_df_L3.loc[j]=Repeat_raw_L3
    return Repeat_df_L3
        #repeatAlign.count_microsatellite(j,V_pattern,window_size,b)

def sliding_window_mini(window_size,V_start,V_length,b,alignAS3,gapAS3):
    Repeat_df_L3 = pd.DataFrame(columns=["Motifs", "Period","RepeatS","Size","Start","End","Position_O","AlignmentS","GapS"])
    for j in range(V_start-window_size+1,V_start+V_length):
        V_pattern = b[j:j+window_size]
        Repeat_raw_L3 = repeatAlign.count_minisatellite(j,V_pattern,window_size,b,V_start,alignAS3,gapAS3)
        if Repeat_raw_L3 != None:
            Repeat_df_L3.loc[j]=Repeat_raw_L3
    return Repeat_df_L3

def search_repeats(df,ref_file,min_size0,max_size0,exact_size0,alignAS0,gapAS0):
    Repeat_df_L0S = pd.DataFrame(columns=["Chromosome","Position","REF","ALT","Genotype","Motifs","Period","Size","Start","End","Repeat_Score","Alignment_Score","Gap_Score","Repeat_GC%"])
    
    # count porprotion of work finished
    #no=0
    #totalno=df.shape[0]
    #process_list = lambda x:[round(y*totalno) for y in x]
    #s = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    #s2 = process_list(s)
    #no_count = 1
    
    position_list = []
    chr_name=""
    for i in range(0,df.shape[0]):
        chromosome =  df.loc[i,"Chromosome"]
        position = df.loc[i,"Position"]
        REF = df.loc[i,"REF"]
        ALT = df.loc[i,"ALT"]
        Genotype = df.loc[i,"Genotype"]

        # skip the searching if current variant exist in the last variants repeat regions, only skip the certain window sizes of last variant repeat region
        # set this part to reduce computational burden,so it will search this variants is exist in previous regions or not
        #if position_list != []:
        #   for previous in position_list:
        #       if int(position) >= previous[2] and int(position) <= previous[3]:
        #           print 1,previous                   
             
        position_list = []

        # read reference
        if chr_name!=chromosome:
            chr_name=chromosome
            Ref=read_reference.read_reference_sequence(ref_file,chr_name)

        # Here, I consider SNV and indel separately, because, SNV, the position in VCF is variant itself, but indel is the reference, 1bp head of variant itself
        # Indel
        if len(REF) != len(ALT):
            Repeat_raw_L0 = catching_variant_Indel(position, REF, ALT, Ref,min_size0,max_size0,exact_size0,alignAS0,gapAS0) # deletion or insertion considered in loop 1
        # SNV or MNV
        else:
            Repeat_raw_L0 = catching_variant_SNV(position, REF, ALT, Ref,min_size0,max_size0,exact_size0,alignAS0,gapAS0)

        Repeat_df_L0 = pd.DataFrame(columns=["Chromosome","Position","REF","ALT","Genotype","Motifs","Period","Size","Start","End","Repeat_Score","Alignment_Score","Gap_Score","Repeat_GC%"])
        if Repeat_raw_L0.empty == False:
            for j in range(0,Repeat_raw_L0.shape[0]):
                if j==0:
                    GC_Percent = Local_GC.GC_content(Ref[Repeat_raw_L0.iloc[j]["Start"].astype(np.int64)-1:Repeat_raw_L0.iloc[j]["End"].astype(np.int64)])
                    item=[chromosome,position,REF,ALT,Genotype,str(Repeat_raw_L0.iloc[j]["Motifs"]),Repeat_raw_L0.iloc[j]["Period"].astype(np.int64),Repeat_raw_L0.iloc[j]["Size"].astype(np.int64),Repeat_raw_L0.iloc[j]["Start"].astype(np.int64),Repeat_raw_L0.iloc[j]["End"].astype(np.int64),round(Repeat_raw_L0.iloc[j]["RepeatS"].astype(np.float64),7),round(Repeat_raw_L0.iloc[j]  ["AlignmentS"].astype(np.float64),7),round(Repeat_raw_L0.iloc[j]["GapS"].astype(np.float64),7),GC_Percent]
                    Repeat_df_L0.loc[j]=item
                    position_list_item = [i,Repeat_raw_L0.iloc[j]["Size"].astype(np.int64),Repeat_raw_L0.iloc[j]["Start"].astype(np.int64),Repeat_raw_L0.iloc[j]["End"].astype(np.int64)]           
                    position_list.append(position_list_item)
                else:
                    GC_Percent = Local_GC.GC_content(Ref[Repeat_raw_L0.iloc[j]["Start"].astype(np.int64)-1:Repeat_raw_L0.iloc[j]["End"].astype(np.int64)])
                    item=["\t","\t","\t","\t","\t",str(Repeat_raw_L0.iloc[j]["Motifs"]),Repeat_raw_L0.iloc[j]["Period"].astype(np.int64),Repeat_raw_L0.iloc[j]["Size"].astype   (np.int64),Repeat_raw_L0.iloc[j]["Start"].astype(np.int64),Repeat_raw_L0.iloc[j]["End"].astype(np.int64),round(Repeat_raw_L0.iloc[j]["RepeatS"].astype(np.float64),7),round(Repeat_raw_L0.iloc[j]   ["AlignmentS"].astype(np.float64),7),round(Repeat_raw_L0.iloc[j]["GapS"].astype(np.float64),7),GC_Percent]
                    Repeat_df_L0.loc[j]=item
                    position_list_item = [i,Repeat_raw_L0.iloc[j]["Size"].astype(np.int64),Repeat_raw_L0.iloc[j]["Start"].astype(np.int64),Repeat_raw_L0.iloc[j]["End"].astype(np.int64)]           
                    position_list.append(position_list_item)
            Repeat_df_L0S = Repeat_df_L0S.append(Repeat_df_L0,ignore_index = True)

            
        else:
            continue
        
        #no=no+1
        #for num in s2:
        #   if no == num:
        #       print str(no_count*10)+"% variant repeat analysis are done"
        #-      no_count = no_count +1

    #print Repeat_df_L0S

    return Repeat_df_L0S



































