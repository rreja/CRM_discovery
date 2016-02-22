import sys, os, re, difflib, math
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np
from operator import itemgetter



def process_file_6(output):
    
    TOMTOM_command = os.path.join(os.path.dirname(os.path.realpath(__file__)),"meme_4.10.2/bin/tomtom")
    # Read the directory containing different states
    for direc in os.listdir(output):
        direc_path = os.path.join(output,direc)
        if not os.path.isdir(direc_path):
            continue
        # directory containing fasta files:
        fasta_dir = os.path.join(direc_path,"peaks")
        outfile = os.path.join(fasta_dir,"local_db.meme")
        out = open(outfile,"w")
        out.write("MEME version 4.4\n\n")
        out.write("ALPHABET= ACGT\n\n")
        out.write("strands: + -\n\n")
        out.write("Background letter frequencies (from uniform background):\n")
        out.write("A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n")
        for fdir in os.listdir(fasta_dir):
            
            if not os.path.isdir(os.path.join(fasta_dir,fdir)):
                continue
            label = fdir.split("_")[0]
            in0 = open(os.path.join(os.path.join(fasta_dir,fdir),"meme.txt"),"rt")
            flag = 0
            flag2 = 0
            motif_count = 1
            for line in in0:
                if line.startswith("log-odds"):
                    flag = 1
                    out.write("MOTIF "+label+"-"+str(motif_count)+"\n\n")
                    out.write(line)
                    motif_count = motif_count + 1
                    continue
                if flag == 1:
                    if line.startswith("-"):
                        flag = 0
                        out.write("\n")
                        continue
                    out.write(line)
                if line.startswith("letter-probability"):
                    flag2 = 1
                    out.write(line)
                    continue
                
                if flag2 == 1:
                    if line.startswith("-"):
                        flag2 = 0
                        out.write("\n")
                        continue
                    out.write(line)
        
        out.close()
        # Run TOMTOM for each meme input on the database
        for fdir in os.listdir(fasta_dir):
            if not os.path.isdir(os.path.join(fasta_dir,fdir)):
                continue
            MEME_FILE = os.path.join(os.path.join(fasta_dir,fdir),"meme.txt")
            os.system("'"+TOMTOM_command+"'"+" -oc "+"'"+os.path.join(fasta_dir,fdir)+"'"+" -min-overlap 5 -thresh 0.05 -dist pearson -evalue -no-ssc "+"'"+MEME_FILE+"'"+" "+"'"+outfile+"'")
    
    
    
    # Create a matrix for each state and cluster to identify motif cluster and motif with a good e-value.
    for direc in os.listdir(output):
        direc_path = os.path.join(output,direc)
        if not os.path.isdir(direc_path):
            continue
        print "Processing state "+direc
        # Create a file that contains the enriched motifs for each state and their E-values
        motif_file = os.path.join(os.path.join(output,direc),"enriched_motifs.txt")
        network_file = os.path.join(os.path.join(output,direc),"network.txt")
        outm = open(motif_file,"w")
        # directory containing fasta files:
        fasta_dir = os.path.join(direc_path,"peaks")
        list1 = []
        evals = {}
        for fdir in os.listdir(fasta_dir):
            if not os.path.isdir(os.path.join(fasta_dir,fdir)):
                continue
            TOMTOM_FILE = os.path.join(os.path.join(fasta_dir,fdir),"tomtom.txt")
            MEME_FILE = os.path.join(os.path.join(fasta_dir,fdir),"meme.txt")
            in0 = open(TOMTOM_FILE,"rt")
            label = fdir.split("_")[0]
            evals = get_evals(MEME_FILE,label,evals)
            for line in in0:
                if line.startswith("#"):
                    continue
                cols = line.rstrip().split("\t")
                query = label+"-"+cols[0]
                ans = cols[1]
                dist = float(cols[3])
                list1.append((query,ans))
                
        
        merged = []
        seen = []
        new_dataset = {}
        for k1,v1 in find_cliques(list1,network_file).items():
            merged = v1
            # if the list belonging to k1 is already been merged then skip it here.
            if k1 in seen:
                continue
            for k2,v2 in find_cliques(list1,network_file).items():
                # if the list belonging to k2 is already been merged then skip it here.
                if k2 in seen:
                    continue
                # Look for overlap between lists and then merger them
                if lists_overlap(merged,v2):
                    merged = list(set(merged)|set(v2))
                    seen.append(k2)
            
                
            new_dataset[k1] =  merged
        
        # Write down the motif for each state with its evalue.
        for k,v in new_dataset.items():
            test_dict = {}
            if len(v) == 1:
                outm.write(v[0]+"\t"+str(evals[v[0]])+"\n")
            else:
                for vals in v:
                    test_dict[vals] = evals[vals]
                sorted_x = sorted(test_dict.iteritems(), key=itemgetter(1),reverse=False)
                outm.write(sorted_x[0][0]+"\t"+str(sorted_x[0][1])+"\n")
        outm.close()
                
        
            
        
                
            
def find_cliques(list1,network_file):
    out = open(network_file,"w")
    d = defaultdict(list)
    for k, v in list1:
        d[k].append(v)
    
    # Printing the network file here, optional code.
    for k,v in d.items():
        for vals in v:
            out.write(k+"\t"+vals+"\n")
    
    
    return d        
    
    
def lists_overlap(a, b):
    return bool(set(a) & set(b))


def get_evals(infile,label,evals):
    in0 = open(infile,"rt")
    motif_count = 1
    for line in in0:
        if line.startswith("log-odds"):
            evals[label+"-"+str(motif_count)] = float(line.rstrip().split("E= ")[1])
            motif_count = motif_count + 1
    
    return(evals)






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
    
    
    (options, args) = parser.parse_args()
    process_file_6(options.output)
    
    
    
if __name__ == "__main__":
    run() 