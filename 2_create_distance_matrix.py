import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
import numpy as np


def process_files(idxdir,options):
    count = 1
    idxData = {}
    for fname in os.listdir(idxdir):
        if fname.endswith("idx") or fname.endswith("tab"):
            input = open(os.path.join(idxdir,fname),"r")
            for line in input:
                if line.startswith("chrom") or line.startswith("#"):
                    continue;
                chrom, start,ftag,rtag,ttag = line.rstrip().split("\t")
                idxData[str(count)+"\t"+chrom+"\t"+start] = int(ttag)
            count = count+1
    print "Done reading all files in memory"
    en_input = open(options.enrcfile,"r")
    for line in en_input:
        
    
            


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
    parser.add_option('-w', action='store', type='int', dest='window',default=200,
                      help='Size of the moving window.Default=200')
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