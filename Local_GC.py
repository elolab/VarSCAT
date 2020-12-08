from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

def GC_content(Sequence):
    GC_value = str(round(GC(Sequence),2))+"%"
    return GC_value
