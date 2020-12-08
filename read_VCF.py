import sys
import pandas as pd
import vcf
import UPS_system
import var_context 
import mut_ref_var
import Complex_variants
import Repeat_context

def Redundant_variant(P,R,A):
    # cutting from rigth 
    while R[-1]==A[-1] and len(R)!=1 and len(A)!=1:
        R=R[:len(R)-1]
        A=A[:len(A)-1]
    # cutting from left
    if len(R) !=1 and len(A) !=1:
        while R[0]==A[0] and len(R)!=1 and len(A)!=1:
            R=R[1:]
            A=A[1:]
            P=str(int(P)+1)
    return P,R,A

def read_vcf(vcf_file,location):
    if len(location) == 0:
        vcf_reader = vcf.Reader(open(vcf_file, "r"))
    else:   
        chromosome_read = location.split(':')[0]
        start_read = location.split(':')[1].split('-')[0]
        end_read = location.split(':')[1].split('-')[1]
        try:
            vcf_reader = vcf.Reader(open(vcf_file, "r")).fetch(str(chromosome_read), int(start_read), int(end_read))
        except ValueError: 
            print ("Error:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")             
            
    sample_name = vcf_reader.samples[0]

    # split multi-alts (how to deal with it remains problem)
    variant_list = []
    
    #count_convenient=1 # convenient use only
    for i in vcf_reader:
        chromosome=str(i.CHROM)
        position = str(i.POS)
        REF = str(i.REF)
        ALT = str(i.ALT[0])
        genotype = str(i.genotype(sample_name)['GT'])
        
        # split multi-alts into single record
        if "," in ALT:
            ALT_multi = ALT.split(",")
            for j in ALT_multi:
                ALT_2 = j
                variant_info = [chromosome,position, REF, ALT_2, genotype]
                variant_list.append(variant_info)
        else:
            variant_info = [chromosome,position, REF, ALT, genotype]
            variant_list.append(variant_info)
        
        # convenient use only
        #count_convenient+=1
        #if count_convenient==50:
        #break

    # Normalize variants: only clean the redundant variants
    # senairo 1: MNP,redundant alleles
    # senairo 2: DEL,redundant alleles
    # senairo 3: INS,redundant alleles
    # redundant alleles:on the right side, there are common letters between REF and ALT
    #                   on the left side, there are common letters between REF and ALT
    variant_list_df = []
    for i in range(0,len(variant_list)):
        if len(variant_list[i][2])!=1 and len(variant_list[i][3])!=1:
            Normalized_variant = Redundant_variant(variant_list[i][1],variant_list[i][2],variant_list[i][3])
            A=Normalized_variant[0]
            B=Normalized_variant[1]
            C=Normalized_variant[2]
            #if len(B)!=1 and len(C)!=1: check if all the variants are normalized
            #    print variant_list[i]   empty output
            item = {"Chromosome":variant_list[i][0],"Position":A,"REF":B,"ALT":C,"Genotype":variant_list[i][4]}
            variant_list_df.append(item)

        else:
            item = {"Chromosome":variant_list[i][0],"Position":variant_list[i][1],"REF":variant_list[i][2],"ALT":variant_list[i][3],"Genotype":variant_list[i][4]}
            variant_list_df.append(item)

    variant_list_df = pd.DataFrame(variant_list_df)
    variant_list_df = variant_list_df[["Chromosome","Position","REF","ALT","Genotype"]]

    return variant_list_df

def complex_module(df,ref,fla,loc,Mut,com,tap,tal,table_def):
    # UPS as a column
    item_UPS = UPS_system.UPS(df,ref)
    df["UPS_L"] = item_UPS[0]
    df["UPS_R"] = item_UPS[1]
    df["UPS_LEFT"] = item_UPS[2]
    df["UPS_RIGHT"] = item_UPS[3]

    # reference seq and mutation seq as input
    # default as 0
    if fla=="1":
        item_varcon = var_context.mutation_seq(df,ref)
        df["ref_sequence"] = item_varcon[0]
        df["mut_sequence"] = item_varcon[1]
    
    if loc!="" and Mut=="1":
        item_mutref = mut_ref_var.mut_ref(df,loc,ref,com,tap,tal,table_def)
    
    # Complex variants as column
    item_complex = Complex_variants.combine_complex(df)
    df["Complex_Cluster"] = item_complex[0]
    df["distance_3'_nearest_Var(bp)"] = item_complex[1]
    del df["UPS_L"]
    del df["UPS_R"]
    
    return df

def repeat_module(df,ref,min_size,max_size,exact_size,alignAS,gapAS):
    # Repeat regions as a column
    item_repeat = Repeat_context.search_repeats(df,ref,min_size,max_size,exact_size,alignAS,gapAS)
    
    return item_repeat


 






























