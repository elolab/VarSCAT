from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Align
from Bio import motifs


def alignment_pattern(PS,SZ,RF,PR,DIR): #PS=pattern start,SZ=search size,RF=ref,PR=potential_R,DIR=direction
    # Setup aligner, values are default from EBI EMBOSS Needle (https://www.ebi.ac.uk/Tools/psa/emboss_needle/)
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5.0
    aligner.mismatch_score = -4.0
    aligner.open_gap_score = -50.0 # do not do gap in this stage

    # we search the possible in a reasonable regions, which are max one pattern size of gap
    list_loc = []
    list_score = []
    list_seq = []
    for i in range(0,SZ):
        if DIR == 1:
            # to right
            potential_R_N = Seq(RF[PS+SZ+i:PS+SZ+SZ+i].upper())
        elif DIR == -1:
            # to left
            potential_R_N = Seq(RF[PS-SZ-i:PS-i].upper())
        score_R = aligner.score(PR, potential_R_N)
        # I remove score balanced with length, I want to use the score with length bias later, to calculate another repeat score which take distance score into account
        # however, I still need to use score balanced with length to condition each alignment (next step), but here I only want to keep raw length bias alognment score for future use # score_R = score_R/len(PR)
        list_loc.append(i)
        list_score.append(score_R)
        list_seq.append(potential_R_N)

    MS_score = max(list_score)
    MS_index = list_score.index(MS_score)
    # distance score should be calcuate here, open gap is -10, extend gap is -1
    MS_distance = list_loc[MS_index]
    if MS_distance==0:
        MS_distance_score = float(MS_distance)
    else:
        MS_distance_score = float(10+(MS_distance-1)*1)
    MS_seq = list_seq[MS_index]

    return MS_score, MS_distance, MS_seq, MS_distance_score



def count_minisatellite(pattern_start, potential_R, search_size, ref, variant_start,alignAS,gapAS):
    # check right pattern
    repeat_time = 1
    repeat_start = pattern_start
    repeat_end = pattern_start + search_size -1

    repeat_start1 = repeat_start
    repeat_end1 = repeat_end

    # for left moving use
    pattern_start1 = pattern_start

    potential_R = Seq(potential_R.upper())
    potential_R0 = potential_R # use to calculate the first
    memory_pattern_list_R = [] # to record every hits that pass the alignment, later to use to generate consensus
    memory_pattern_list_R.append(potential_R)

    # to remember the indels of alignment, later use this to see how consensus the repeat regions are
    memory_distance_list = []
    # to remember the scores of alignment (which based on alignment setting, it will reflect mismatch and indel level)
    # later use this to see how consensus the repeat regions are
    memory_score_list = []

    # to right
    Align_inner = alignment_pattern(pattern_start,search_size,ref,potential_R,1)  # 1: here use as direction, 1 means right, -1 mean left
    # here, we still calcuate score balanced with length to condition each alignment
    # I took 4bp as example, I think MMMW is ok
    while float(Align_inner[0])/len(potential_R) >=2.75:  #Threshold:MMMW, (5x3-4)/4, Align_inner[0] is the most similar alignment score
        if "N" not in Align_inner[2]:
            memory_pattern_list_R.append(Align_inner[2]) # Align_inner[2] is the pattern that matched
            memory_distance_list.append(Align_inner[3]) # Align_inner[3] is the distance score
            memory_score_list.append(Align_inner[0])
            if len(memory_pattern_list_R) > 4: # more than 4 nucleotides ATGC
                m = motifs.create(memory_pattern_list_R)
                potential_R = m.consensus
            else:
                potential_R = potential_R
            pattern_start = pattern_start+search_size + Align_inner[1] # Align_inner[1] is the actually gap between each pattern
            repeat_time = repeat_time +1
            repeat_end = repeat_end + search_size + Align_inner[1]
            Align_inner = alignment_pattern(pattern_start,search_size,ref,potential_R,1)
        else:
            break

    # to left
    Align_inner = alignment_pattern(pattern_start1,search_size,ref,potential_R,-1)  # 1: here use as direction, 1 means right, -1 mean left
    while float(Align_inner[0])/len(potential_R) >=2.75:
        if "N" not in Align_inner[2]:
            memory_pattern_list_R.append(Align_inner[2])
            memory_distance_list.append(Align_inner[3])
            memory_score_list.append(Align_inner[0])
            if len(memory_pattern_list_R) > 4:
                m = motifs.create(memory_pattern_list_R)
                potential_R = m.consensus
            else:
                potential_R = potential_R
            pattern_start1 = pattern_start1 - search_size -  Align_inner[1]
            repeat_time = repeat_time +1
            repeat_start = repeat_start -search_size -  Align_inner[1]
            Align_inner = alignment_pattern(pattern_start1,search_size,ref,potential_R,-1)
        else:
            break

    if repeat_start != repeat_start1 or repeat_end != repeat_end1:
        # according repeat defination
        if repeat_time >=5:
            # For the concept of repeat, the first four alignment block should be the variant block align to consensus pattern
            aligner = Align.PairwiseAligner()
            aligner.mode = "global"
            aligner.match_score = 5.0
            aligner.mismatch_score = -4.0
            aligner.open_gap_score = -50
            score_R = aligner.score(potential_R0, potential_R)
            memory_score_list.append(score_R)
            for i in range(1,4):
                score_R = aligner.score(memory_pattern_list_R[i], potential_R)
                memory_score_list.append(score_R)
            del memory_score_list[1:4]
            # average distance between each repeat unit (bp)
            distance_score = float(sum(memory_distance_list))/len(memory_distance_list)
            # average alignment score of each repeat unit
            alignment_score = float(sum(memory_score_list))/len(memory_score_list)
            # average alignment score of each base
            alignment_score1 = alignment_score/len(potential_R)
            #alignment_score2 = float(5-alignment_score1)/9 # proportion of mismatches of one pattern
            # here, I think distance can be gaps of alignment, which is -1 for aligment score. Divided by units length is going to get a univeral parameters
            #repeat_score = float(alignment_score-distance_score)/len(potential_R)
            #12.3 new socre, I think the whole repeat region as a whole
            #sum(memory_score_list) is all the match=1, mismatch=0 in this regions, sum(memory_distance_list) is all the open gap = -1, extend gap=-0.5 in this region
            #then sum them all together is the regions alignment for the pattern, this_score/len(regions) ----> length free alignment score
            # Or minimize the influence of length
            #repeat_score = float(sum(memory_score_list)-sum(memory_distance_list))/float(len(potential_R))#float(repeat_end-repeat_start+1)
            # the two version repeat score is not equal, the former one is larger
            repeat_score = (float(sum(memory_score_list)-sum(memory_distance_list))/(repeat_end-repeat_start+1))*repeat_time
            # theoretical the minimum score: MMMWGGGMMMWGGGMMMWGGGMMMWGGGMMMW M:match, W:mismatch, D:gap
            # 19.11, here I changed the criteria
            if repeat_score >= alignAS*repeat_time and distance_score <= gapAS: # 3 is the minimum per base quality as default, 5 is average more than 0.5 gap
                # repeat_start is the python corrdinates start, repeat_end is the python coordinates end
                # all +1 because this is python position, for function "alignment_pattern",[start:end] is not the same as repeat_start and repeat_end
                # repeat_end I -1 in function "count_minisatellite", due to I want to use the last base of repeat region as end.
                # but in "alignment_pattern", I didn't -1. So in the end, repeat_start and repeat_end need to +1
                return potential_R.upper(),repeat_time,repeat_score,len(potential_R.upper()),repeat_start+1,repeat_end+1,variant_start,alignment_score1,distance_score


def count_microsatellite(pattern_start, potential_R, search_size, ref, variant_start):
    # check right pattern
    repeat_time = 1
    repeat_start = pattern_start
    repeat_end = pattern_start + search_size-1

    repeat_start1 = repeat_start
    repeat_end1 = repeat_end

    pattern_start1 = pattern_start

    potential_R = Seq(potential_R.upper())
    # to right
    while potential_R == Seq(ref[pattern_start+search_size:pattern_start+search_size+search_size].upper()):
        pattern_start = pattern_start+search_size
        repeat_time = repeat_time +1
        repeat_end = repeat_end + search_size
    # to left
    while potential_R == Seq(ref[pattern_start1-search_size:pattern_start1].upper()):
        pattern_start1 = pattern_start1 - search_size
        repeat_time = repeat_time +1
        repeat_start = repeat_start -search_size
    if repeat_start != repeat_start1 or repeat_end != repeat_end1:
        # according repeat defination
        if repeat_time >=5:
            # print potential_R, repeat_time, repeat_start, repeat_end
            # python 0-based, convert back to VCF 1-based, so start+1. Python do not include last index, so repeat_end -1 is the last based on repeat for 0-based, for 1-based should +1, -1+1=0
            # print potential_R.upper(), repeat_time, repeat_start+1, repeat_end, "microsatellite"
            # fill up the column with minisatellite, although empty
            distance_score = 0
            alignment_score = 5
            repeat_score = repeat_time*5 # if simply the repeatS in minisatellite, the formula looks like this
            return potential_R.upper(),repeat_time,repeat_score,len(potential_R.upper()),repeat_start+1,repeat_end+1,variant_start,alignment_score,distance_score












































