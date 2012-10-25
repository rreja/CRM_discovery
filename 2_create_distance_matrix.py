import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
import numpy as np
from random import randint


def process_files(idxdir,options):
    count = 1
    idxData = {}
    encData = {}
    allData = {}
    filehash = {}
    taglist = []
    enriched = open(options.enrcfile,"r")
    for line in enriched:
        chrom,junk,junk,start,end,junk,junk,junk,junk  = line.rstrip().split("\t")
        encData[chrom+":"+start+":"+end] = 1
    
    # Reading all idx files in the directory one by one
    for fname in os.listdir(idxdir):
        if fname.endswith("idx") or fname.endswith("tab"):
            input = open(os.path.join(idxdir,fname),"r")
            for line in input:
                if line.startswith("chrom") or line.startswith("#"):
                    continue;
                chrom, start,ftag,rtag,ttag = line.rstrip().split("\t")
                idxData[chrom+"\t"+start] = int(ttag)
                
        # Calculating the tags present in the enrcihed regions.
        for key, val in encData.items():
            l = key.split(":")
            for i in range(int(l[1]),int(l[2])):         
                idx = l[0]+"\t"+str(i)
                if idx in idxData:
                    taglist.append(idxData[idx])
                else:
                    taglist.append(0)
                    
            allData[str(count)+":"+key] = taglist
            taglist = []
        filehash[count] = fname
        count = count+1
    
    create_matrix_for_regions(allData)               

def create_matrix_for_regions(allData):
    for majorkey, majorval in allData.items():
        for minorkey, minorval in allData.items():
            compute_distance(majorkey,minorkey)
            
                        
def print_dict(dictionary):
    for key, val in dictionary.items():
        print key,val


usage = '''
input_paths may be:
- a single file.

example usages:
python 2_create_distance_matrix.py [options] /dir_to_idx_files/
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-e', action='store',  dest='enrcfile',
                      help='File containing the enriched regions.')
    parser.add_option('-w', action='store', type='int', dest='window',default=100,
                      help='Size of the moving window.Default=100')
    parser.add_option('-s', action='store', type='int', dest='stp_size',default=1,
                      help='Slide the window with given bp length in an interval.Defalut=1')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    if not os.path.isdir(args[0]):
        print "Input data direcotry containing all idx files."
    else:
        process_files(args[0],options)
    
if __name__ == "__main__":
    run() 