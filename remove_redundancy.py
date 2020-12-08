from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Align
from Bio import motifs

## df.size-- total number of cell, df.ndim--number of dimensions
## This part should contain two part:1. merge motifs with same size, 2. kick out redundancy

def clean_df(df): # double clean repeats, use pattern simularity
    df = df.sort_values('RepeatS',ascending=False) # Based on which criteria, it remain discussion 3017430  769454
    df = df.reset_index()
    drop_list = []
    for i in range(0,df.shape[0]-1): # df.shape[0]--row number, df.shape[1]--column number
        for j in range(i+1,df.shape[0]):
            region_ij = float(max(df.iloc[i]['End'],df.iloc[j]['End'])-min(df.iloc[i]['Start'],df.iloc[j]['Start'])+1)
            # overlap (concept of union) of two repeat regions is less than 0.3
            # redundant_score = float((min(df.iloc[i]['End'],df.iloc[j]['End'])-max(df.iloc[i]['Start'],df.iloc[j]['Start'])+1)*2)/region_ij
            redundant_score = float(min(df.iloc[i]['End'],df.iloc[j]['End'])-max(df.iloc[i]['Start'],df.iloc[j]['Start'])+1)/region_ij
            if redundant_score > 0.3:
                drop_list.append(j)

    drop_list = set(drop_list)
    drop_list = list(drop_list)
    df = df.drop(df.index[drop_list])
    return df

# drop a row: df.drop(['Cochice', 'Pima'])
# drop a column: df.drop('reports', axis=1)
# Drop a row by row number (in this case, row 3): df.drop(df.index[2])
# find duplicated line: df[df.duplicated(['Size'],keep=False)]
