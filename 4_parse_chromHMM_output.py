import sys, os, re, difflib, fnmatch
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np
import pybedtools
from operator import itemgetter

def process_file(options):
    
    print "Creating directory for each state."
    # Get the file names to process
    for fname in os.listdir(options.output):
        if not (fname.endswith(".bed") or fname.endswith(".txt")):
            continue
        if fname.startswith("emissions"):
            emissions_fname =  fname
        if fname.endswith("_segments.bed"):
            segment_fname = fname
    
    in0 = open(os.path.join(options.output,segment_fname),"rt")
    segments = defaultdict(list)
    for line in in0:
        cols = line.rstrip().split("\t")
        segments[cols[3]].append(line.rstrip())
    
    in0.close()
    # Create a direcotry for each state and write a segment file for each state.
    for k,v in segments.items():
        
        # Create new directory
        new_dir = os.path.join(options.output,k)
        if not os.path.exists(new_dir): os.makedirs(new_dir)
        
        # write segment file.
        out = open(os.path.join(new_dir,k+"_segments.gff"),"w")
        for vals in v:
            cols = vals.split("\t")
            start = int(cols[1]) + 1
            end = int(cols[2])
            out.write(cols[0]+"\t.\t.\t"+str(start)+"\t"+str(end)+"\t.\t.\t.\t.\n")
        out.close()
    
    
    print "Getting enriched factors in each state!"
    # Write a file containing the factors enriched in a particular state.
    in1 = open(os.path.join(options.output,emissions_fname),"rt")
    for line in in1:
        cols = line.rstrip().split("\t")
        if line.startswith("state"):
            label = cols[1:]
            continue
        state = cols[0]
        pvals = [float(x) for x in cols[1:]]
        outdir = os.path.join(options.output,"E"+state)
        
        # Create a peak direcotry in each folder
        peak_dir = os.path.join(outdir,"peaks")
        if not os.path.exists(peak_dir): os.makedirs(peak_dir)
        
        out = open(os.path.join(outdir,"factors.txt"),"w")
        for i,j in zip(label,pvals):
            if j >= options.filter:
                
                out.write(i+"\n")
                # Writing the peak files that intersect with enriched segments.
                infile = return_filename(options.peakdir,i+"*")
                sloped = return_filename(outdir,"*_segments.gff")
                inter = os.path.join(peak_dir,infile)
                pybedtools.BedTool(os.path.join(options.peakdir,infile)).intersect(os.path.join(outdir,sloped),u=True).saveas(inter)
                
                # Get top 500 peak-pairs from each and slop distance upstream/downstream.
                sorted_inter = os.path.join(peak_dir,os.path.splitext(infile)[0]+"_top500.gff")
                get_top500(inter,sorted_inter,options)
                os.system("rm "+"'"+inter+"'")
                
                # extract fasta sequences
                fastaseq = os.path.join(peak_dir,os.path.splitext(infile)[0]+"_sequences.fa")
                write_out = open(fastaseq,"w")
                a = pybedtools.BedTool(sorted_inter)
                
                for line in open(a.sequence(fi=options.fasta).seqfn).readlines():
                    write_out.write(line)
                write_out.close()
                
                # Remove un-necessary files.
                os.system("rm "+"'"+sorted_inter+"'")
                
                
        out.close()
    

def get_top500(infile,outfile,options):
    in0 = open(infile,"rt")
    out = open(outfile,"w")
    data = {}
    for line in in0:
        if line.startswith("#"):
            out.write(line)
            continue
        cols = line.rstrip().split("\t")
        data[line] = float(cols[5])
    
    in0.close()
    sorted_x = sorted(data.iteritems(), key=itemgetter(1),reverse=True)
    count = 1
    for j in sorted_x:
        if count > 500:
            break
        cols = j[0].rstrip().split("\t")
        start = int(cols[3])-options.dist
        end = int(cols[4])+options.dist
        if start <=0 :
            continue
        out.write(cols[0]+"\t"+cols[1]+"\t"+cols[2]+"\t"+str(start)+"\t"+str(end)+"\t"+cols[5]+"\t"+cols[6]+"\t"+cols[7]+"\t"+cols[8]+"\n")
        count = count + 1
    
    out.close()
    
def return_filename(outdir,pattern):
    for file in os.listdir(outdir):
        if fnmatch.fnmatch(file, pattern):
            return(file)
   
                
        
   




usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python 4_parse_chromHMM_output.py  [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-o', action='store', type='string', dest='output',
                      help='ChromHMM Output directory')
    parser.add_option('-f', action='store', type='string', dest='filter',default = 0.18,
                      help='Threshold for emission probability to consider enrichment in a state, default = 0.18')
    parser.add_option('-p', action='store', type='string', dest='peakdir', 
                      help='Directory containing original peak call files.')
    parser.add_option('-s', action='store', type='int', dest='dist',default = 30, 
                      help='Upstream/downstream distance from mid-point.')
    parser.add_option('-i', action='store', type='string', dest='fasta', 
                      help='Reference FASTA file.')
    
    
    (options, args) = parser.parse_args()
    
    ## output dir for ChromHMM output
    #hmmdir = os.path.join(options.dir,"_chromHMM_OUTPUT")
    #if not os.path.exists(hmmdir): os.makedirs(hmmdir)
    
    process_file(options)
    
    
    
if __name__ == "__main__":
    run() 