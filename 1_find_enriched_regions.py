import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
import numpy as np
import pybedtools

def process_files(dataDir,options):
    input = open(options.reffile,'rt')
    # Create a concatenation of all peak-pair gffs to calculate an overall pp file.
    outfile2 =  os.path.join(dataDir,"all.pp.gff")
    #alldata  = os.path.join(dataDir,"*.gff")
    #os.system("cat "+alldata+" >"+outfile2)
    mline = ""
    #outfile1 =  os.path.join(dataDir,"tmp.genome.gff") 
    outfile3 = os.path.join(dataDir,"all.enriched.gff")
    outfile4 = os.path.join(dataDir,"all.background.gff")
    #tmp_genome = open(outfile1,'r+')
    for line in input:
        if line.startswith("#"):
            continue
        chrom,end = line.split("\t")
        #get_coordinates(chrom,1,int(end),tmp_genome,options)
        mline = mline+get_coordinates(chrom,1,int(end),options)
        
    # Creating bedtool object for all 1000 bp genome coordinates.
    all_genome = pybedtools.BedTool(mline,from_string=True)
    # Applying the intersection function to get the enriched and non-enriched regions.
    

    all_genome.intersect(outfile2,u=True).saveas(outfile3) # Gives all the unique chromosomal 1000bp coordinate that are enriched.
    all_genome.intersect(outfile2,v=True).saveas(outfile4) # Gives all the 1000 bp coordinate that are background.
    os.system("rm "+outfile2)
    
    

def get_coordinates(chrom,start,end,options):
    line = ""
    for i in range(start,end,options.step_size):
       
        if not i+999 > end:
            #genome_writer.write(chrom+"\t.\t.\t"+str(i)+"\t"+str(i+999)+"\t.\t.\t.\t.\n")
            line = line+chrom+"\t.\t.\t"+str(i)+"\t"+str(i+999)+"\t.\t.\t.\t.\n"
            
        else:
            #genome_writer.write(chrom+"\t.\t.\t"+str(i)+"\t"+str(end)+"\t.\t.\t.\t.\n")
            line = line+chrom+"\t.\t.\t"+str(i)+"\t"+str(end)+"\t.\t.\t.\t.\n"
    return(line)
        
    




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
    parser.add_option('-s', action='store', type='int', dest='step_size',default=1000,
                      help='Step size to use to fragment genome.')
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