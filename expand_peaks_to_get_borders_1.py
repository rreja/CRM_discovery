import sys, os, math, operator
from optparse import OptionParser , IndentedHelpFormatter
from itertools import tee, izip, islice
from scipy import stats
import numpy as np
import pybedtools
from operator import itemgetter
import pysam
 



def process_file_1(options,outdir):
    
    genome_length = get_start_end(options.gfile)
    # Calculation of total mappable windows in the genome
    total_windows =  options.gsize/options.window
    #print "Total number of mappable windows in the genome = "+str(total_windows)
    
    BAM_files = {}
    # Read through the BAM files
    for fname in os.listdir(options.BAMdir):
        if not fname.endswith(".bam"):
            continue
        bam_name = fname.rstrip().split("_")[0]
        BAM_files[bam_name] = os.path.join(options.BAMdir,fname)
    
    
    # Read through the peak-pair directory.
    for fname in os.listdir(options.peakDir):
        if not fname.endswith(".gff"):
            continue
        
        # Get the factors name [MAKE SURE YOUR FILENAME STARTS WITH FACTORNAME FOLLOWED BY "_" ]
        outfile = os.path.join(outdir,os.path.splitext(fname)[0]+"_expanded.gff")
        label = fname.split("_")[0]
        print "Processing "+fname
        in0 = open(os.path.join(options.peakDir,fname),"rt")
        peak_info = {}
        for line in in0:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            peak_info[cols[8]] = cols[0]+":"+cols[3]+":"+cols[4]
        
        in0.close()
        
        # Send the peak information and the corresponding BAM file to the function
        if label not in BAM_files:
            print "Please check the naming schema of your BAM file for factor ="+fname
            print "Both the BAM file and the peak file name should start with the factor name followed by underscore '_' "
            sys.exit(1)
            
        new_peak_borders = expand_peak_borders(peak_info,BAM_files[label],options,genome_length,total_windows)
        
        #print "Writing the peaks info for "+fname
        out = open(outfile,"w")
        for k,v in new_peak_borders.items():
            vals = v.split(":")
            out.write(vals[0]+"\t.\t"+k+"\t"+vals[1]+"\t"+vals[2]+"\t.\t.\t.\t"+k+"\n")
        out.close()
        
            
        
def expand_peak_borders(peak_info,bamFile,options,genome_length,total_windows):
    
    samfile = pysam.AlignmentFile(bamFile, "rb")
    # Calculation of background read distribution
    total_read = get_total_reads_from_BAM(samfile)
    expected_reads_per_window = float(total_read)/total_windows
    if expected_reads_per_window < 1:
        expected_reads_per_window = 1
    #print "Total reads in file = "+str(total_read)
    #print "Expected reads per window = "+str(expected_reads_per_window)
    p = 0.5
    
    new_peak_borders = {}
    for k,v in peak_info.items():
        chrom = v.split(":")[0]
        start = int(v.split(":")[1])
        
        # Go upstream from start.
        flag = 1
        new_start = start - options.window
        end = start
        while(flag):
            
            # Check if the start cooridnate does not become negative
            if not check_start_end(chrom,new_start,end,genome_length):
                upstream_border = 1
                break
            
            # Get tags in the interval
            sum_tags = get_sum_in_interval(samfile,chrom,new_start,end)
            n = sum_tags+expected_reads_per_window
            pval = stats.binom.sf(sum_tags,n,p)
            if pval > options.pval:
                upstream_border = new_start
                break
            
            # Create new coordinates
            old_start = new_start
            new_start = old_start - options.window
            end = old_start
        
        # Go downstream from start.  
        flag = 1
        new_end = start + options.window
        while(flag):
            
             # Check if the start cooridnate does not become negative
            if not check_start_end(chrom,start,new_end,genome_length):
                downstream_border = genome_length[chrom]
                break
            
            # Get tags in the interval
            sum_tags = get_sum_in_interval(samfile,chrom,start,new_end)
            n = sum_tags+expected_reads_per_window
            pval = stats.binom.sf(sum_tags,n,p)
            if pval > options.pval:
                downstream_border = new_end
                break
            
            # Create new coordinates
            old_end = new_end
            new_end = old_end + options.window
            start = old_end
            
        new_peak_borders[k] = chrom+":"+str(upstream_border)+":"+str(downstream_border)
        
    return(new_peak_borders)

def get_sum_in_interval(samfile,chrom,start,end):
    sum_tags = 0
    for read in samfile.fetch(chrom,start,end):
        
        # Check if the read is paired or single end:
        if read.is_paired:
            
            # Only consider read1:
            if not read.is_read1:
                continue
            # Only consider proper-pair:
            if not read.is_proper_pair:
                continue
            # Find out if the read is mapped to the sense strand or the reverse to get the 5p end
            if read.is_reverse:
                read_5p = read.reference_end
            else:
                # Need to convert the BED 5p coordinate to gff as my intervals from peak file are in gff format
                read_5p = read.reference_start+1    
            # Check if the 5p end lies in your defined interval
            if (read_5p < start or read_5p >= end):
                continue
            sum_tags = sum_tags+1
            
        else:
            if read.is_unmapped:
               continue
             # Find out if the read is mapped to the sense strand or the reverse to get the 5p end
            if read.is_reverse:
                read_5p = read.reference_end
            else:
                # Need to convert the BED 5p coordinate to gff as my intervals from peak file are in gff format
                read_5p = read.reference_start+1
            # Check if the 5p end lies in your defined interval
            if (read_5p < start or read_5p >= end):
                continue
            sum_tags = sum_tags+1
            
    
    return(sum_tags)
    
def get_total_reads_from_BAM(samfile):
    total_reads = 0
    for read in samfile.fetch():
        if read.is_paired:
            if not read.is_proper_pair:
                continue
            if not read.is_read1:
                continue
            total_reads = total_reads + 1
        else:
            if read.is_unmapped:
               continue
            total_reads = total_reads + 1
    
    return(total_reads)
        

def check_start_end(chrom,start,end,genome_length):
    if start > 0 and end <= genome_length[chrom]:
        return(True)
    else:
        return(False)
    
def get_start_end(gfile):
    chrom_ends = {}
    in0 = open(gfile,"rt")
    for line in in0:
        cols = line.rstrip().split("\t")
        chrom_ends[cols[0]] = int(cols[1])
    return(chrom_ends)

def windows(iterable, size):
    iters = tee(iterable, size)
    for i in xrange(1, size):
        for each in iters[i:]:
            next(each, None)
    return izip(*iters)



usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python 1_expand_peaks_to_get_borders.py [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-p', action='store', type='string', dest='peakDir',
                      help='The directory containing all the peaks/peak-pairs calls.')
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
        
    
    process_file_1(options,outdir)
    
    
if __name__ == "__main__":
    run() 