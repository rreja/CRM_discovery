import sys, os, math, operator
from optparse import OptionParser , IndentedHelpFormatter
from itertools import tee, izip, islice
from scipy import stats
import numpy as np
import pybedtools
from operator import itemgetter
import pysam
from expand_peaks_to_get_borders_1 import process_file_1
from create_and_run_chomHMM_2 import process_file_2



def process_file(options):
    
    outdir = os.path.join(options.peakDir,"_expanded_regions")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    ## Pipeline step-1: Expand peak-pairs to get border
    print "STEP-1: Expanding peaks to get border."
    #process_file_1(options,outdir)
    print "STEP-1 completed!"
    
    
    # Pipeline step-2: Create chromHMM input and run chromHMM
    ### Input directory for chromHMM input
    chromhmm_outdir = os.path.join(outdir,"_chromHMM_INPUT")
    if not os.path.exists(chromhmm_outdir): os.makedirs(chromhmm_outdir)
    
    ### Output directory for ChromHMM output
    chromhmm_indir = os.path.join(outdir,"_chromHMM_OUTPUT")
    if not os.path.exists(chromhmm_indir): os.makedirs(chromhmm_indir)
    
    print "STEP-2: Creating chromHMM Input."
    process_file_2(options,options.peakDir,options.gfile,chromhmm_outdir,chromhmm_indir)
    print "STEP-2 completed!"
    
    # Pipeline step-3: Create chromHMM input and run chromHMM


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python Run_pipeline.py [OPTIONS]
**** Remember that your BAM file and your peak file should start with "factorname_"
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
                      help='Window size, default = 20')
    parser.add_option('-g', action='store', type='string', dest='gfile',
                      help='File containing the chromosome number and length.')
    parser.add_option('-s', action='store', type='int', dest='gsize',default=11332237,
                      help='Mappable genome size: sg11 = 11,332,237 (default), hg19=248,988,565, mm9=2,178,433,024 for read length = 36. Refer PMID:22276185')
    parser.add_option('-v', action='store', type='float', dest='pval',default=0.05,
                      help='P-value cutoff for significant enrichment over background, default = 0.05')
    parser.add_option('-l', action='store', type='int', dest='ilen', default = 200,
                      help='Length of chromHMM segmentation interval, default = 200')
    parser.add_option('-n', action='store', type='int', dest='state',default = 12, 
                      help='Number of states, default = 12')
    
    (options, args) = parser.parse_args()
    process_file(options)
    
    
if __name__ == "__main__":
    run() 