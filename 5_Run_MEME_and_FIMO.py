import sys, os, re, difflib, fnmatch
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np
import pybedtools
from operator import itemgetter

def process_file(options):
    
    MEME_command = os.path.join(os.path.dirname(os.path.realpath(__file__)),"meme_4.10.2/bin/meme")
    FIMO_command = "'"+os.path.join(os.path.dirname(os.path.realpath(__file__)),"meme_4.10.2/bin/fimo")+"'"
    
    # Read the directory containing the fasta files
    for direc in os.listdir(options.output):
        direc_path = os.path.join(options.output,direc)
        if not os.path.isdir(direc_path):
            continue
        # directory containing fasta files:
        fasta_dir = os.path.join(direc_path,"peaks")
        for fname in os.listdir(fasta_dir):
            if not fname.endswith(".fa"):
                continue
            fasta_file = os.path.join(fasta_dir,fname)
            outdir = os.path.join(fasta_dir,os.path.splitext(fname)[0])
            
            # Run MEME
            #print "'"+MEME_command+"'"+" -dna -mod zoops -nmotifs "+str(options.nmotifs)+" -minsites "+str(options.minsites)+" -revcomp -oc "+"'"+outdir+"'"+" "+fasta_file
            
            os.system("'"+MEME_command+"'"+" -dna -mod zoops -nmotifs "+str(options.nmotifs)+" -minsites "+str(options.minsites)+" -revcomp -oc "+"'"+outdir+"'"+" "+"'"+fasta_file+"'")
            
            
            # Run FIMO
            #os.system(FIMO_command+" --parse-genomic-coord  --thresh "+str(options.thresh)+" --motif 1 -oc "+outdirFIMO+" "+os.path.join(outdir,"meme.html")+" "+options.inFIMO)
        sys.exit(1)
    
    
    
    



                
        
   




usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python 5_Run_MEME_and_FIMO.py  [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-o', action='store', type='string', dest='output',
                      help='ChromHMM Output directory')
    parser.add_option('-f', action='store', type='int', dest='nmotifs',default = 3,
                      help='Number of motifs, defualt = 3')
    parser.add_option('-p', action='store', type='int', dest='minsites',default = 50, 
                      help='Minimum sites to contain motif, default = 50.')
    parser.add_option('-s', action='store', type='float', dest='thresh',default = 0.001, 
                      help='Threshold for FIMO motif detection.')
    parser.add_option('-i', action='store', type='string', dest='fasta', 
                      help='Reference FASTA file.')
    
    
    (options, args) = parser.parse_args()
    process_file(options)
    
    
    
if __name__ == "__main__":
    run() 