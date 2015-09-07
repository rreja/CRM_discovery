import sys, os, pybedtools, random
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
import numpy as np
from operator import itemgetter
import math

def process_file(options):
    
    states = {}
    in0 = open(options.file,"rt")
    total_intervals = 0
    for line in in0:
        cols = line.rstrip().split("\t")
        no_of_intervals = (int(cols[2]) - int(cols[1]))/200
        total_intervals = total_intervals + no_of_intervals
        if cols[3] in states:
            states[cols[3]] = states[cols[3]] + no_of_intervals
        else:
            states[cols[3]] =  no_of_intervals
        
    
    in0.close()
    
    print total_intervals
    bins_cal_after = 0
    for k,v in states.items():
        print "No of bins = "+str(v)
        print "Percentage of bins in state "+k+" "+str((float(v)/total_intervals)*100)
        bins_cal_after = bins_cal_after + v
    print bins_cal_after


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python Parse_hmm_output_for_matrix.py
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-i', action='store', type='string', dest='file',
                      help='File containing segments information in BED format.')
    #parser.add_option('-j', action='store', type='string', dest='tabDir',
    #                  help='tab/idx directory.')
    #parser.add_option('-r', action='store', type='string', dest='ref',
    #                  help='+1 nucleosome coordinates')
    
    
    
    (options, args) = parser.parse_args()       
    
    #outfile = os.path.join(os.path.dirname(options.hmm),"Mean_vector_matrix.txt")    
    process_file(options)
    
    
    

    
if __name__ == "__main__":
    run() 