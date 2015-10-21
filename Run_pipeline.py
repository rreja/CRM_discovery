import sys, os, math, operator
from optparse import OptionParser , IndentedHelpFormatter
from itertools import tee, izip, islice
from scipy import stats
import numpy as np
import pybedtools
from operator import itemgetter
import pysam
from expand_peaks_to_get_borders_1 import *



def process_file(options,outdir):
    
    process_file_1(options,outdir)
    print "done"
    
 


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python Run_pipeline.py [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-p', action='store', type='string', dest='peakDir',
                      help='The directory containing the peaks call files in gff format.')
    parser.add_option('-i', action='store', type='string', dest='BAMdir',
                      help='The directory containing the BAM files.')
    parser.add_option('-w', action='store', type='int', dest='window',default=20,
                      help='Window size')
    parser.add_option('-g', action='store', type='string', dest='gfile',
                      help='File containing the chromosome number and length.')
    parser.add_option('-s', action='store', type='int', dest='gsize',default=11332237,
                      help='Mappable genome size: sg11 = 11,332,237 (default), hg19=248,988,565, mm9=2,178,433,024 for read length = 36. Refer PMID:22276185')
    parser.add_option('-v', action='store', type='float', dest='pval',default=0.05,
                      help='P-value cutoff for significant enrichment over background, default = 0.05')
    
    (options, args) = parser.parse_args()
    
    outdir = os.path.join(options.peakDir,"_expanded_regions")
    if not os.path.exists(outdir): os.makedirs(outdir)
        
    
    process_file(options,outdir)
    
    
if __name__ == "__main__":
    run() 