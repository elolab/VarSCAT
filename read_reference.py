# read reference
import sys
from Bio import SeqIO

def read_reference_sequence(input_file,chromosome):
    try:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            if seq_record.id==chromosome:
                sequence_chr=str(seq_record.seq)
        return sequence_chr
    except UnboundLocalError:
        print ("Error:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
