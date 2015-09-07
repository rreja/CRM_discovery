import sys, os, pybedtools, random
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
from math import asinh, log
import math



def process_file(options,outfile):
    in0 = open(options.ref,"rt")
    reference = {}
    for line in in0:
        if line.startswith("#") or line.startswith("chrom"):
            continue
        cols = line.rstrip().split("\t")
        id_key = cols[0]+"-"+cols[3]+"-"+cols[4]
        #reference[cols[0]+":"+cols[8]] = int((int(cols[3])+int(cols[4]))/2)
        reference[cols[0]+":"+id_key] = int((int(cols[3])+int(cols[4]))/2)
    
    in0.close()
    
    std_data = defaultdict(list)
    
    header = "GENES"
    for fname in os.listdir(options.tabDir):
        if not fname.endswith(".tab"):
            continue
        in1 = open(os.path.join(options.tabDir,fname),"rt")
        print "Processing "+fname
        header = header+"\t"+fname.split("_")[0]
        idxData = {}
        for line in in1:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            idxData[cols[0]+":"+cols[1]] = float(cols[2]) + float(cols[3])
        
        in1.close()
        
        std_data = calculate_std(options,reference,std_data,idxData)
    

    
    
    
    out = open(outfile,"w")
    out.write(header+"\n")
    #out.write(header+"\trandom\n")
    for k,v in std_data.items():
        line = k
        for vals in v:
            line = line+"\t"+str(vals)
        out.write(line+"\n")
    
    out.close()
        

def calculate_std(options,reference,std_data,idxData):
    for k,v in reference.items():
        chrom = k.split(":")[0]
        gene = k.split(":")[1]
        start = v - options.dist
        end = v + options.dist
        variance = []
        for j in range(start,end+1):
            key = chrom+":"+str(j)
            if key not in idxData:
                continue
            variance.append(idxData[key]*(j-v)*(j-v))
        
        std_data[gene].append(math.pow(sum(variance)/len(variance),0.5))
        
    
    return(std_data)
            
            
        
            
            
    



usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python slop_and_intersectBed.py /usr/local/tab_files
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-r', action='store', type='string', dest='ref',
                      help='Reference files in gff format.')
    parser.add_option('-i', action='store', type='string', dest='tabDir',
                      help='Direcotry containing tab files.')
    parser.add_option('-d', action='store', type='int', dest='dist',default = 500,
                      help='Upstream/downstream distance from the reference. Default = 500')
    

    (options, args) = parser.parse_args()
    outfile = os.path.join(options.tabDir,"deviations_from_ref.txt")
    process_file(options,outfile)
    
    
    
    
if __name__ == "__main__":
    run() 