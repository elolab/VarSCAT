def combine_complex(df):    
    left_coordinates = []
    right_coordinates = []
    for i in range(0,df.shape[0]):
        left_coordinates.append(int(df.loc[i,"UPS_L"]))
        right_coordinates.append(int(df.loc[i,"UPS_R"]))

    left_coordinates = left_coordinates[1:]
    right_coordinates = right_coordinates[:len(right_coordinates)-1]
    
    k=0
    m=0
    n=0
    e=0
    f=0
    CLS=[]
    distance_to_3T = []
    for i in range(0,len(left_coordinates)):
        #ambiguous_value_A = abs(len(df.loc[i,"REF"])-len(df.loc[i,"ALT"]))
        #ambiguous_value_B = int(df.loc[i,"UPS_R"])-int(df.loc[i,"UPS_L"])
        # Indels and the indel regions is not variant length -> ambiguous
        #ambiguous_value = (ambiguous_value_A!=ambiguous_value_B)
        distance_to_3T.append(left_coordinates[i]-right_coordinates[i])
        if left_coordinates[i]<=(right_coordinates[i]+1): #and ambiguous_value == True:
            # control clsuter number
            if m!=n:
                k=k+1
                CLS1 = df.loc[i,"Chromosome"]+"_"+df.loc[i,"Position"]+"_"+df.loc[i,"REF"]+"_"+df.loc[i,"ALT"]+"_"+df.loc[i,"Genotype"]+"_"+"CLS"+str(k)
                CLS2 = df.loc[i+1,"Chromosome"]+"_"+df.loc[i+1,"Position"]+"_"+df.loc[i+1,"REF"]+"_"+df.loc[i+1,"ALT"]+"_"+df.loc[i,"Genotype"]+"_"+"CLS"+str(k)
            if CLS1 not in CLS:
                CLS.append(CLS1)
            if CLS2 not in CLS:    
                CLS.append(CLS2)
                f=f+1
                
            m=n
        else:
            if e!=f:
                e=f
                continue
            else:
                CLS.append("")
                e=f
            n=n+1
    CLS.append("")
    distance_to_3T.append("")
    
    return CLS,distance_to_3T




