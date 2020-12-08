import sys
import read_reference

#Here, I need to build a UPS system of all the variants, in order to find the variants range, to build the cluster

def cutting_right(start_C,pattern_C,pattern_D,b): # move variant to left most position
    if len(pattern_C)>len(pattern_D):
        variants_C = pattern_C.upper()[1:] # only get variant itself, not include REF
    elif len(pattern_C)<len(pattern_D):
        variants_C = pattern_D.upper()[1:]
    suffix_coordinate = len(variants_C) # last letter of variants
    suffix = variants_C[suffix_coordinate-1]
    flank_left_coordinate = int(start_C)-1   # start_C should be the position of variant in VCF, which is the reference coordinates of REF. In O-based python, it should subtract 1
    # flank_left = b[flank_left_coordinate].upper() this part is removed, in ordre to deal with complex variants such as TCG A--
    flank_left = pattern_D.upper()[0].upper()
    variants_pattern=variants_C

    if b[int(start_C)-1].upper() == pattern_D.upper()[0].upper():
        while suffix == flank_left:
            variants_pattern = variants_pattern[len(variants_pattern)-1]+variants_pattern[:len(variants_pattern)-1]
            suffix =variants_pattern[suffix_coordinate].upper()
            flank_left_coordinate = flank_left_coordinate-1
            flank_left = b[flank_left_coordinate].upper()
        
        if len(pattern_C)>len(pattern_D):
            # change python 0-based into reference 1-based, this is REF reference coordinate in VCF files; +1 again, because I want to use the variant itself
            N_start_C = flank_left_coordinate+1+1 
        elif len(pattern_C)<len(pattern_D):
            N_start_C = flank_left_coordinate+1
    else:
        if len(pattern_C)>len(pattern_D):
            # change python 0-based into reference 1-based, this is REF reference coordinate in VCF files; +1 again, because I want to use the variant itself
            N_start_C = flank_left_coordinate+1+1 
        elif len(pattern_C)<len(pattern_D):
            N_start_C = flank_left_coordinate+1

    return N_start_C


def cutting_left(start_C,pattern_C,pattern_D,b): # move variant to right most position
    if len(pattern_C)>len(pattern_D):
        variants_C = pattern_C.upper()[1:] # only get variant itself, not include REF
    elif len(pattern_C)<len(pattern_D):
        variants_C = pattern_D.upper()[1:]
    prefix_coordinate = 0
    prefix = variants_C[prefix_coordinate]
    # the right flank for deletion and insertion are different corrdinates. CGAT C--T(T=position+lenCGA-1 is the flank right). C--T CGAT (T=position+1-1, T is the flank right)
    if len(pattern_C)>len(pattern_D):
        flank_right_coordinate = int(start_C)+len(pattern_C)-1  
    elif len(pattern_C)<len(pattern_D):
        flank_right_coordinate = int(start_C)
    flank_right = b[flank_right_coordinate].upper() 
    variants_pattern=variants_C
    
    if b[int(start_C)-1].upper() == pattern_D.upper()[0].upper():
        while prefix == flank_right:
            variants_pattern = variants_pattern[1:]+variants_pattern[0]
            #print variants_pattern
            prefix =variants_pattern[prefix_coordinate].upper()
            flank_right_coordinate = flank_right_coordinate+1
            flank_right = b[flank_right_coordinate].upper()

        N_end_C = flank_right_coordinate # change python 0-based into reference 1-based, this is the last letter of variant in reference coordinate
    else:
        N_end_C = flank_right_coordinate
    
    return N_end_C

def UPS(df,ref_file):
    LEFT=[]
    RIGHT=[]
    LEFT_2=[]
    RIGHT_2=[]
    chr_name=""
    for i in range(0,df.shape[0]):
        # read reference
        chromosome=df.loc[i,"Chromosome"]
        if chr_name!=chromosome:
            chr_name=chromosome
            ref=read_reference.read_reference_sequence(ref_file,chr_name)

        if len(df.loc[i,"REF"])>len(df.loc[i,"ALT"]):
            variants1 = cutting_right(df.loc[i,"Position"],df.loc[i,"REF"],df.loc[i,"ALT"],ref)
            LEFT.append(variants1)
            variants2 = cutting_left(df.loc[i,"Position"],df.loc[i,"REF"],df.loc[i,"ALT"],ref)
            RIGHT.append(variants2)
            if variants1 != int(df.loc[i,"Position"])+1 or variants2 != int(df.loc[i,"Position"])+len(df.loc[i,"REF"])-1:
                LEFT_2.append(variants1)
                RIGHT_2.append(variants2)
            else:
                LEFT_2.append("")
                RIGHT_2.append("")

        elif len(df.loc[i,"REF"])<len(df.loc[i,"ALT"]):
            variants1 = cutting_right(df.loc[i,"Position"],df.loc[i,"REF"],df.loc[i,"ALT"],ref)
            LEFT.append(variants1)
            variants2 = cutting_left(df.loc[i,"Position"],df.loc[i,"REF"],df.loc[i,"ALT"],ref)
            RIGHT.append(variants2)
            if variants1 != int(df.loc[i,"Position"]) or variants2 != int(df.loc[i,"Position"]):
                LEFT_2.append(variants1)
                RIGHT_2.append(variants2)
            else:
                LEFT_2.append("")
                RIGHT_2.append("")

        elif len(df.loc[i,"REF"])==len(df.loc[i,"ALT"]):
            LEFT.append(int(df.loc[i,"Position"]))
            RIGHT.append(int(df.loc[i,"Position"])+len(df.loc[i,"REF"])-1)
            LEFT_2.append("")
            RIGHT_2.append("")

    return LEFT,RIGHT,LEFT_2,RIGHT_2





















