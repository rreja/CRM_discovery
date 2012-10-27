import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
import pybedtools

def process_files(dataDir,options):
    chr_length = {}
    mline = ""
    input = open(options.reffile,'rt')
    for line in input:
        chrom,maxlen = line.rstrip().split("\t")
        chr_length[chrom] = maxlen
    
    # Create a newDir output, so that we don't concatenate the results into query:
    newDir = os.path.join(dataDir,"output")
    if not os.path.exists(newDir):
            os.mkdir(newDir)
            
    # Create a concatenation of all peak-pair gffs to calculate an overall pp file.
    outfile2 =  os.path.join(newDir,"all.pp.gff")
    alldata  = os.path.join(dataDir,"*.gff")
    os.system("cat "+alldata+" >"+outfile2)
    
    # Creating outfiles for downstream analysis.
    outfile3 = os.path.join(newDir,"all.enriched.gff")
    outfile4 = os.path.join(newDir,"all.background.gff") 
    #outfile100 = os.path.join(newDir,"all.merged100.gff")
    
    ## Trying to merge overlapping segments in  gff file.The output is a bed file. So I need to add +1 to the start.
    b = pybedtools.BedTool(outfile2).merge(d=options.dist)
    #pybedtools.BedTool(outfile2).merge(d=options.dist).saveas(outfile100)
      
    for i in b:  
        if abs(int(i[2]) - (int(i[1])+1)) < options.window:
            # Extend short regions.
            start =  (int(i[1])+1) - options.window/2
            end   =  int(i[2]) + options.window/2
            if start <= 0 or end > int(chr_length[i[0]]):
                if start <=0:
                    mline = mline+i[0]+"\t.\t.\t"+str(1)+"\t"+str(end)+"\t.\t.\t.\t.\n"
                else:
                    mline = mline+i[0]+"\t.\t.\t"+str(start)+"\t"+str(chr_length[i[0]])+"\t.\t.\t.\t.\n"
            else:
                mline = mline+i[0]+"\t.\t.\t"+str(start)+"\t"+str(end)+"\t.\t.\t.\t.\n"
            
        else:
            start = int(i[1])+1
            end   = int(i[2])
            if start <=0:
                mline = mline+i[0]+"\t.\t.\t"+str(1)+"\t"+str(end)+"\t.\t.\t.\t.\n"
            elif end > int(chr_length[i[0]]):
                mline = mline+i[0]+"\t.\t.\t"+str(start)+"\t"+str(chr_length[i[0]])+"\t.\t.\t.\t.\n"
            else:
                mline = mline+i[0]+"\t.\t.\t"+str(start)+"\t"+str(end)+"\t.\t.\t.\t.\n"
            
    pybedtools.BedTool(mline,from_string=True).saveas(outfile3)
    # complement will be in BED format.
    pybedtools.BedTool(mline,from_string=True).complement(g=options.reffile).saveas(outfile4)
    os.system("rm "+outfile2)
    
   

#def get_coordinates(chrom,start,end,options):
#    line = ""
#    for i in range(start,end,options.step_size):
#       
#        if not i+999 > end:
#            #genome_writer.write(chrom+"\t.\t.\t"+str(i)+"\t"+str(i+999)+"\t.\t.\t.\t.\n")
#            line = line+chrom+"\t.\t.\t"+str(i)+"\t"+str(i+999)+"\t.\t.\t.\t.\n"
#            
#        else:
#            #genome_writer.write(chrom+"\t.\t.\t"+str(i)+"\t"+str(end)+"\t.\t.\t.\t.\n")
#            line = line+chrom+"\t.\t.\t"+str(i)+"\t"+str(end)+"\t.\t.\t.\t.\n"
#    return(line)
#        
    




usage = '''
input_paths may be:
- a single file.

example usages:
python 1_find_enriched_regions.py [options] /dir_to_data_files/
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-g', action='store',  dest='reffile',
                      help='File containing chromosome length')
    # Here the window size should always be less then or equal to the "d". Since two regions will always be "d" distance apart,
    # if we extend regions by window/2 bp (for regions less then the scanning window size), we can still have non-overlapping regions.
    parser.add_option('-w', action='store', type='int', dest='window',default=100,
                      help='Size of scanning window to be used in downstream analysis. Default=100.')
    parser.add_option('-d', action='store', type='int', dest='dist',default=100,
                      help='Merge contigs "dist" bp apart into a single regions. Default=100')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    try:
        import pybedtools
    except ImportError:
        print "You need to install pybedtools."
        sys.exit(1)
     
    if not os.path.isdir(args[0]):
        print "Input data direcotry containing all gff files."
    else:
        process_files(args[0],options)
    
if __name__ == "__main__":
    run() 