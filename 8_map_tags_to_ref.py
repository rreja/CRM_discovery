import sys, os, re, difflib, math
from optparse import OptionParser , IndentedHelpFormatter
from collections import defaultdict
from operator import itemgetter
import pybedtools
import pysam



def process_file(options):
    
    bams = {}
    # Creating list of BAM files    
    for fname in os.listdir(options.bam):
        if not fname.endswith(".bam"):
            continue
        labs = fname.split("_")[0]
        bams[labs] = os.path.join(options.bam,fname)
        
    # Read the directory containing different states
    for direc in os.listdir(options.output):
        direc_path = os.path.join(options.output,direc)
        if not os.path.isdir(direc_path):
            continue
        
        
        
        # Get information about the motif
        Tss_file = os.path.join(direc_path,"enriched_TSS.gff")
        
        for fname in os.listdir(os.path.join(direc_path,"_candidate_refs")):
            if not fname.endswith(".gff"):
                continue
            label = fname.split("-")[0]
            # Create directory for motifs
            
            #Create dicrecotories for CDT files.
            outdir = os.path.join(direc_path,"_CDT_"+label)
            if not os.path.exists(outdir): os.makedirs(outdir)
            
            ref_file = os.path.join(os.path.join(direc_path,"_candidate_refs"),fname)
            bam_file = bams[label]
            samfile = pysam.AlignmentFile(bam_file, "rb")
            cdt_file = {}
            in0 = open(ref_file,"rt")
            
            for line in in0:
                if line.startswith("chrom") or line.startswith("#"):
                    continue
                cols = line.rstrip().split("\t")
                chrom = cols[0]
                
                if cols[6] == "+":
                    start = int(cols[3]) - options.up
                    end = int(cols[3]) + options.down
                else:
                    start = int(cols[4]) - options.down
                    end = int(cols[4]) + options.up
                
                key = cols[0]+":"+str(start)+":"+str(end)+":"+cols[6]
                cdt_file[key] = get_tags_in_interval(samfile,chrom,start,end)
                
            in0.close()
            write_sense_antisense_file(cdt_file,outdir,label,options)
    
            
        
def write_sense_antisense_file(cdt_file,outdir,label,options):
    out_sense = open(os.path.join(outdir,label+"_sense.cdt"),"w")
    out_anti = open(os.path.join(outdir,label+"_antisense.cdt"),"w")
    write_header(out_sense,options)
    write_header(out_anti,options)
    for k,v in cdt_file.items():
        chrom = k.split(":")[0]
        start = int(k.split(":")[1])
        end = int(k.split(":")[2])
        strand = k.split(":")[3]
        sense = []
        anti = []
        
        # Create sense and anti lists to store the cooridnates at each position
        for j in range(start,end+1):
            if strand == "+":
                key_pos = chrom+":"+str(j)+":+"
                key_neg = chrom+":"+str(j)+":-"
                if key_pos in v:
                    sense.append(v[key_pos])
                else:
                    sense.append(0)
                
                if key_neg in v:
                    anti.append(v[key_neg])
                else:
                    anti.append(0)
                    
            elif strand == "-":
                key_pos = chrom+":"+str(j)+":+"
                key_neg = chrom+":"+str(j)+":-"
                if key_neg in v:
                    sense.append(v[key_neg])
                else:
                    sense.append(0)
                
                if key_pos in v:
                    anti.append(v[key_pos])
                else:
                    anti.append(0)
                    
        # Write the sense list as it is but for TSS on antisense strand flip the coordinates.
        if strand == "+":
            line = k+"\t1"
            for j in sense:
                line = line+"\t"+str(j)
            out_sense.write(line+"\n")
            line = k+"\t1"
            
            for k in anti:
                line = line+"\t"+str(k)
            out_anti.write(line+"\n")
        
        elif strand == "-":
            line = k+"\t1"
            for j in reversed(sense):
                line = line+"\t"+str(j)
            out_sense.write(line+"\n")
            
            line = k+"\t1"
            for k in reversed(anti):
                line = line+"\t"+str(k)
            out_anti.write(line+"\n")
            
    out_sense.close()
    out_anti.close()


def get_tags_in_interval(samfile,chrom,start,end):
    data = {}
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
                strand = "-"
               
                
            else:
                # Need to convert the BED 5p coordinate to gff as my intervals from peak file are in gff format
                read_5p = read.reference_start+1
                strand = "+"
               
                
            # Check if the 5p end lies in your defined interval
            if (read_5p < start or read_5p >= end):
                continue
            
            
            key = chrom+":"+str(read_5p)+":"+strand
            if key in data:
                data[key] = data[key] + 1
            else:
                data[key] = 1
            
            
        else:
            if read.is_unmapped:
               continue
             # Find out if the read is mapped to the sense strand or the reverse to get the 5p end
            if read.is_reverse:
                read_5p = read.reference_end
                strand = "-"
            else:
                # Need to convert the BED 5p coordinate to gff as my intervals from peak file are in gff format
                read_5p = read.reference_start+1
                strand = "+"
            # Check if the 5p end lies in your defined interval
            if (read_5p < start or read_5p >= end):
                continue
            
            key = chrom+":"+str(read_5p)+":"+strand
            if key in data:
                data[key] = data[key] + 1
            else:
                data[key] = 1

    return(data)
        
        
        
            
def  write_header(out,options):
    line = "Uniqe ID\tGWEIGHT"
    for j in range(-1*(options.up),options.down+1):
        line = line+"\t"+str(j)
    out.write(line+"\n")
    
    
    
    


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python 8_maps_tags_to_ref.py  [OPTIONS]
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-o', action='store', type='string', dest='output',
                      help='ChromHMM output directory')
    parser.add_option('-i', action='store', type='string', dest='bam', 
                      help='Directory containing BAM files.')
    parser.add_option('-u', action='store', type='int', dest='up',default = 500,
                      help='Upstream distance from ref = 500.')
    parser.add_option('-d', action='store', type='int', dest='down',default = 500,
                      help='Downstream distance from ref = 500.')
    
    
    (options, args) = parser.parse_args()
    process_file(options)
    
    
    
if __name__ == "__main__":
    run() 