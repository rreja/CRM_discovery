import sys, os, operator, random, time
from optparse import OptionParser , IndentedHelpFormatter
from itertools import izip, cycle, tee

def process_file(idxdir,options):
    idxData = {}
    filehash = {}
    count = 1
    # Reading all idx files in the directory one by one
    for fname in os.listdir(idxdir):
        if fname.endswith("idx") or fname.endswith("tab"):
            input = open(os.path.join(idxdir,fname),"r")
            for line in input:
                if line.startswith("chrom") or line.startswith("#"):
                    continue;
                chrom, start,ftag,rtag,ttag = line.rstrip().split("\t")
                # Appending the count string in the key to recognize each tag file seperately in the same idxData dict
                idxData[str(count)+":"+chrom+":"+start] = int(ttag)
            filehash[count] = fname
            count = count+1
    
    # Create a separate CRM file for each CRM and then give the dir as input. Needs to be implemented.
    for outFile in os.listdir(options.outDir):
        output_folder = os.path.join(options.outDir,os.path.basename(outFile).split('.')[0]) 
        if not os.path.exists(output_folder): os.makedirs(output_folder)
        for no, fname in filehash.items():
            fwrite = os.path.join(output_folder,os.path.basename(fname).split('.')[0]+".cdt")
            out = open(fwrite,"w")
            input2 = open(os.path.join(options.outDir,outFile),"r")
            for line in input2:
                key,windowNo = line.rstrip().split("\t")
                chrom,start,end = get_actual_cord(key,windowNo,options)
                printLine = key+"\t."
                for j in range(start,end):
                    lookup_key = str(no)+":"+chrom+":"+str(j)
                    if lookup_key in idxData:
                        printLine = printLine+"\t"+str(idxData[lookup_key])
                    else:
                        printLine = printLine+"\t0"
                out.write(printLine+"\n")
                out.flush()
                
                
           




def get_actual_cord(key,windowNo,options):
    tmp = key.split(":")
    start = int(tmp[1]) + options.bins*int(windowNo)
    end = start + options.windowLength
    return(tmp[0],start,end)










usage = '''
input_paths may be:

example usages:
python create_CDT_from_output.py -o <dir_containing_CRMs> <dir_to_tag_files>
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description
    
  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-w', action='store', type='int', dest='windowLength',default=100,
                      help='[OPTIONAL]: Size of the moving window.Should be same as the one used in previous step.Default=100')
    parser.add_option('-b', action='store', type='int', dest='bins',default=5,
                      help='[OPTIONAL]: Sub bin within a window. Should be same as the one used in previous step.Default=5')
    parser.add_option('-o', action='store', dest='outDir',
                      help='[MANDATORY]: Output File containing the motif locations and offsets.')
    (options, args) = parser.parse_args()
       
    if not args:
        parser.print_help()
        sys.exit(1)
        
    process_file(args[0],options)   
    
    
if __name__ == "__main__":
    run() 