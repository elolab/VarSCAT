import sys, getopt
import read_VCF
import read_reference

# read reference
def reference_fetch(reference_file):
    ref_seq=read_reference.read_reference_sequence(reference_file)
    return ref_seq

def main():
    vcf_input=''
    reference_f=''
    location=''
    output_complex='output_complex'
    output_repeat='output_repeat'
    flank="1"
    Mut_seq="0"
    Complement="0"
    Transcript="0"
    Translate="0"
    table_def=""
    rmin=1
    rmax=11 # python not included last coordinate, so actually here will be 10 in range(Rmin,Rmax)
    esize=3
    alignment_score=3
    gap_score=10
    # I need these require inputs: -V,--vcf; -R,--reference; -L, --location; -M, --output_complex; -N, --output-repeat
    # these without inputs: -c,--complex; -r, --repeat; -h
    try:
        # define main input
        # 'flank=' is an user options in complex module, to output the mutation sequence context around variants. If not provided, then it should be default that only output UPS variants
        options, remainder = getopt.gnu_getopt(sys.argv[1:], 'crh',['complex','repeat','help','vcf=', 'reference=','location=','output_complex=','output_repeat=','flank=','mut_seq=','complement=','transcript=','translate=','translate_table=','rmin=','rmax=','esize=','alignment_score=','gap_score='])
        # output main help page
        if len(sys.argv[1:])==1 and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
            print('''Help main:
VarCont v1.0
-c,--complex: complex variants analysis module
-r,--repeat: repeat region variants analysis module
-h,--help: help page

Two modules can be used together, the commom parameters '--vcf','--reference' and '--location' should be only announced once.
Note: the program was built on 0 based coordination system
''')    
            sys.exit(2)
        
        # output complex module help
        if ('-c' in sys.argv[1:] or '--complex' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]):
            print('''Help page of complex variants module:
Required parameters:
--vcf: input VCF file 
--reference: input reference sequencing file

Optional parameters:
--location: genome location need to be parsed (the VCF file should be indexed if -L is activated, tbi file of VCF is required, empty value: annotate all variants in VCF)
--output_complex: complex variants module output prefix
--flank: output the flank bases of variants (default=1, equal to True)
--mut_seq: output the reference and mutated sequence based on variants (default=0,0:off,1:on. Note: Only valid when --location is given)
--complement: output the complement sequence of mutated sequence (default=0, Note: Only valid when --mut_seq is given)
--transcript:output the transcription of mutated sequence (default=0, Note: Only valid when --mut_seq is given)
--translate:output the translation of mutated sequence (default=0, Note: Only valid when --mut_seq is given)
--translate_table: the translation table (default=1, Standard table. More info:https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11. Note: Only valid --translation is given)

--help: help page
''')
            sys.exit(2)

        # output repeat module help
        if ('-r' in sys.argv[1:] or '--repeat' in sys.argv[1:]) and ('-h' in sys.argv[1:] or '--help' in sys.argv[1:]): 
            print('''Help page of repeat region variants module:
Required parameters:
--vcf: input VCF file 
--reference: input reference sequencing file

Optional parameters:
--location: genome location need to be parsed (the VCF file should be indexed if -L is activated, tbi file of VCF is required, empty value: annotate all variants in VCF)
--output_repeat: repeat variants module output prefix

Advanced parameters:
--rmin: the minimun size of repeat pattern (default=1)
--rmax: the maximum size of repeat pattern (default=10, larger size will increase the running time)
--esize: the maximum size of repeat pattern for exact match (default=3, the value should between rmin and rmax)
--alignment_score: the minimum average alignment score of each base in repeat regions (default=3, maximum 5)
--gap_score: the maximum average gap of each repeat unit (default=1)

--help: help page
''')
            sys.exit(2)
        
        # prepare reference & normalized vcf:
        for opt,arg in options:
            if opt in ('--reference'):
                reference_f=arg 
            elif opt in ('--location'):
                location=arg
            elif opt in ('--vcf'):
                vcf_input=arg   
        pre_df = read_VCF.read_vcf(vcf_input,location)
        try:
            pre_df = read_VCF.read_vcf(vcf_input,location)
            print ("Normalized VCF ready")
        except:
            print ("VCF file error!")
            sys.exit(2)

        # complex module input fetch
        if ('-c' in sys.argv[1:] or '--complex' in sys.argv[1:]):
            print ("processing complex variants module")
            for opt,arg in options:
                if opt in ('--output_complex'):
                    output_complex=arg
                elif opt in ('--flank'):
                    flank=arg
                elif opt in ('--mut_seq'):
                    Mut_seq=arg
                elif opt in ('--complement'):
                    Complement=arg
                elif opt in ('--transcript'):
                    Transcript=arg
                elif opt in ('--translate'):
                    Translate=arg   
                elif opt in ('--translate_table'):
                    table_def=arg
                
           try:
                complex_result = read_VCF.complex_module(pre_df,reference_f,flank)
                # write into csv
                output_file_complex = output_complex+'.csv'
                complex_result.to_csv (output_file_complex, index = False, header=True)
                print "Complex variants module finished"
           except:
                print "Complex variants module error!!!!"

        # repeat module input fetch
        if ('-r' in sys.argv[1:] or '--repeat' in sys.argv[1:]):
            print ("processing repeat regions variants module")
            for opt,arg in options:
                if opt in ('--output_repeat'):
                    output_repeat=arg
                elif opt in ('--rmin'):
                    rmin=int(arg)
                elif opt in ('--rmax'):
                    rmax=int(arg)
                elif opt in ('--esize'):
                    esize=int(arg)
                elif opt in ('--alignment_score'):
                    alignment_score=float(arg)
                elif opt in ('--gap_score'):
                    gap_score=float(arg)*10
            try:
                repeat_result = read_VCF.repeat_module(pre_df,reference_f,rmin,rmax,esize,alignment_score,gap_score)
                # write into csv
                output_file_repeat = output_repeat+'.csv'
                repeat_result.to_csv (output_file_repeat, index = False, header=True)
                print ("Repeat regions variants module finished")
            except:
                print ("Repeat regions variants module error!!!!")


    # error raise
    except getopt.GetoptError:
        print('''Getopt error! help:
VarCont v1.0
-c,--complex: complex variants analysis module
-r,--repeat: repeat region variants analysis module
-h,--help: help page

Note: the program was built on 0 based coordination system
''')
        sys.exit(2)



if __name__ == "__main__":
    main()
















































