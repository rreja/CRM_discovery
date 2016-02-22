import sys, os, re, difflib, fnmatch
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np
import pybedtools
from operator import itemgetter

def process_file_3(chromHMM_input,chromHMM_output):
    count_ones = {}
    count_zeros = {}
    
    # Get the file names to process
    for fname in os.listdir(chromHMM_output):
        if not (fname.endswith(".bed") or fname.endswith(".txt")):
            continue
        if fname.startswith("emissions"):
            emissions_fname =  fname
        
    
    in0 = open(os.path.join(chromHMM_output,emissions_fname),"rt")
    ems = defaultdict(list)
    state_order = []
    for line in in0:
        cols = line.rstrip().split("\t")
        if cols[0].startswith("state"):
            header = cols[1:]
            for key in header:
                count_ones[key] = 0
                count_zeros[key] = 0
            continue
        
        for k,v in zip(header,cols[1:]):
            ems[k].append(float(v))
            
            
    
    in0.close()
    
    # Reading the input directory.
    
    for fname in os.listdir(chromHMM_input):
        if not fname.endswith(".txt"):
            continue
        in0 = open(os.path.join(chromHMM_input,fname),"rt")
        for line in in0:
            if line.startswith("BY4741"):
                continue
            if re.match(r"^[A-Za-z]",line):
                header = line.rstrip().split("\t")    
                continue
            cols = line.rstrip().split("\t")
            for k,v in zip(header,cols):
                if int(v) == 1:
                    count_ones[k] = count_ones[k] + 1
                elif int(v) == 0:
                    count_zeros[k] = count_zeros[k] + 1
    
    # Calculating the threshold values for each factor belonging to one state.
    states = defaultdict(list)
    state_vals = defaultdict(list)
    
    for k,v in ems.items():
        pct_genome = float(count_ones[k])/(count_ones[k]+count_zeros[k])
        state = 1
        header = "States"
        for vals in v:
            fold_change = vals/(pct_genome)
            #state_vals[state].append(fold_change)
            
            #if fold_change >= 4:
            if vals >= 0.18:
                states[state].append(k)
                #print vals,state
            header = header+"\t"+str(state)
            state = state + 1

    
    out = open(os.path.join(os.path.dirname(chromHMM_output),"emission_threshold_distribution.txt"),"w")
    out.write(header+"\n")
    for j in range(0,len(ems.keys())):
        line = "lines"
        for k,v in state_vals.items():
            line = line+"\t"+str(v[j])
        out.write(line+"\n")
    
    
    return(states)
    #for k,v in states.items():
    #    print k,v, len(v)
                    
            
            
            
   
    
    
    
    
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
    parser.add_option('-i', action='store', type='string', dest='input',
                      help='ChromHMM Input directory')
    
    
    
    (options, args) = parser.parse_args()
    process_file_3(options.input,options.output)
    
    
    
if __name__ == "__main__":
    run() 