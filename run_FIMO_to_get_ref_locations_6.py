import sys, os, re, difflib, math, fnmatch
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np
from operator import itemgetter
import pybedtools



def process_file(options):
    
    FIMO_command = os.path.join(os.path.dirname(os.path.realpath(__file__)),"meme_4.10.2/bin/fimo")
    # Read the directory containing different states
    for direc in os.listdir(options.output):
        direc_path = os.path.join(options.output,direc)
        if not os.path.isdir(direc_path):
            continue
        
        # Create directory for motifs
        outdir = os.path.join(direc_path,"_candidate_refs")
        if not os.path.exists(outdir): os.makedirs(outdir)
        
        # Get information about the motif
        infile = os.path.join(direc_path,"enriched_motifs.txt")
        motif_info = parse_motif_file(infile)
        sloped = return_filename(direc_path,"*_segments.gff")
        # directory containing fasta files:
        fasta_dir = os.path.join(direc_path,"peaks")
        
        for motif in motif_info:
            motif_name = motif.split("-")[0]
            motif_number = int(motif.split("-")[1])
            for fdir in os.listdir(fasta_dir):
                if not os.path.isdir(os.path.join(fasta_dir,fdir)):
                    continue
                label = fdir.split("_")[0]
                if label != motif_name:
                    continue
                
                MEME_FILE = os.path.join(os.path.join(fasta_dir,fdir),"meme.txt")
                infile = os.path.join(outdir,"fimo.gff")
                inter = os.path.join(outdir,motif+".gff")
                os.system("'"+FIMO_command+"'"+" --parse-genomic-coord  --thresh "+str(options.thresh)+" --motif "+str(motif_number)+" -oc "+"'"+outdir+"'"+" "+"'"+MEME_FILE+"'"+" "+options.fasta)
                pybedtools.BedTool(infile).intersect(os.path.join(direc_path,sloped),u=True).saveas(inter)
                os.system("rm "+"'"+infile+"'")
    
            





def parse_motif_file(infile):
    in0 = open(infile,"rt")
    motifs = []
    for line in in0:
        cols = line.rstrip().split("\t")
        motifs.append(cols[0])
        
    return(motifs)
            
        
def return_filename(outdir,pattern):
    for file in os.listdir(outdir):
        if fnmatch.fnmatch(file, pattern):
            return(file)          
            







usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python 7_run_FIMO_to_get_ref_locations.py  [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-o', action='store', type='string', dest='output',
                      help='ChromHMM Output directory')
    parser.add_option('-i', action='store', type='string', dest='fasta', 
                      help='Reference FASTA file.')
    parser.add_option('-t', action='store', type='float', dest='thresh',default = 0.001,
                      help='threshold for FIMO.default = 0.001.')
    
    
    (options, args) = parser.parse_args()
    process_file(options)
    
    
    
if __name__ == "__main__":
    run() 