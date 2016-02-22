import sys, os, re, difflib
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np


def process_file_2(options,peakDir,sg11,outdir,hmmdir):
    
    order = []
    peak_data = defaultdict(list)
    ## Reading the peak-pair directory
    for fname in os.listdir(peakDir):
        if not fname.endswith(".gff"):
            continue
        label = fname.split("_")[0]
        order.append(label)
        in0 = open(os.path.join(peakDir,fname),"rt")
        for line in in0:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            peak_data[label+":"+cols[0]].append(cols[3]+":"+cols[4])
        
        in0.close()
    
            
    # Read the chromosome length file
    in0 = open(sg11,"rt")
    for line in in0:
        cols = line.rstrip().split("\t")
        chrom = cols[0]
        start = 1
        end = int(cols[1])
        outfile = os.path.join(outdir,chrom+"_binary.txt")
        out = open(outfile,"w")
        out.write("BY4741\t"+chrom+"\n")
        out.write("\t".join(order)+"\n")
        
        print "Processning chromosome = "+chrom
        for j in range(start,end,options.ilen):
            int_start = j
            int_end = j + options.ilen - 1
            if int_end > end:
                continue
            out.write(get_value(order,peak_data,chrom,int_start,int_end)+"\n")
        
        out.close()
            
    print "Completed creating input!"


    ## Run ChromHMM
    # Edit file paths to accomodate blank spaces in file name.
    chromHMM = os.path.join(os.path.dirname(os.path.realpath(__file__)),"ChromHMM")
    chromHMM_command = "'"+chromHMM+"/ChromHMM.jar'"
    outdir_modified = "'"+outdir+"'"
    hmmdir_modified = "'"+hmmdir+"'"
    
    #print "java -mx1600M -jar "+chromHMM_command+" LearnModel "+outdir+" "+hmmdir+" "+str(options.state)+" sg11"
    os.system("java -mx1600M -jar "+chromHMM_command+" LearnModel "+outdir_modified+" "+hmmdir_modified+" "+str(options.state)+" sg11 ")
    
    print "Completed Running chromHMM, output in "+hmmdir
    



def get_value(order,peak_data,chrom,start,end):
    
    binary_values = []
    for factor_name in order:
        value = 0
        key = factor_name+":"+chrom
        if not key in peak_data:
            binary_values.append(str(value))
            continue
        for vals in peak_data[key]:
            peak_start = int(vals.split(":")[0])
            peak_end = int(vals.split(":")[1])
            if end >= peak_start and start <= peak_end:
                value = 1
                break
            #if peak_start >= start and peak_end <= end:
            #    value = 1
            #    break
        
        binary_values.append(str(value))
    
    return("\t".join(binary_values))
            
                
                
        
   




usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python create_chromHMM_input_file.py  [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-r', action='store', type='string', dest='ref',
                      help='Reference file containing the chormosome and length.')
    parser.add_option('-i', action='store', type='string', dest='dir',
                      help='Directory containing the peak-pair files.')
    parser.add_option('-l', action='store', type='int', dest='ilen', default = 200,
                      help='Length of interval')
    parser.add_option('-s', action='store', type='int', dest='state',default = 12, 
                      help='Number of states, default = 10')
    
    
    (options, args) = parser.parse_args()
    # output dir for ChromHM input
    outdir = os.path.join(options.dir,"_chromHMM_INPUT")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    # output dir for ChromHMM output
    hmmdir = os.path.join(options.dir,"_chromHMM_OUTPUT")
    if not os.path.exists(hmmdir): os.makedirs(hmmdir)
    
 
    process_file_2(options,options.dir,options.ref,outdir,hmmdir)
    
    
    
if __name__ == "__main__":
    run() 