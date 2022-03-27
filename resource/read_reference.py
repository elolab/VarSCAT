import sys
from pyfaidx import Fasta

def read_reference_sequence(input_file,chromosome):
	try:
		reference_seq = Fasta(input_file,rebuild=False)
		sequence_chr=reference_seq[chromosome][:].seq
		return sequence_chr
	except KeyError:
		print ("Error2:reference file chromosome name, vcf record chromosome name, location chromosome name should be consistent (either all chrX or only X)")
		sys.exit(2)
